"""Module for querying Ensembl API.

This module also contains code to read the VCF file and prepare
the data for querying the Ensembl API.
"""

import contextlib
import csv
import io
import json
import logging
import pathlib
import time
import typing

import httpx
import pysam

from varona import enum

logger = logging.getLogger("varona.ensembl")


class Assembly(enum.CiStrEnum):
    """Enum for the assembly choice.

    There are only two choices for now, GRCh37 and GRCh38.
    """

    GRCh37 = enum.auto()  # aka hg19 or human_g1k_v37
    GRCh38 = enum.auto()  # aka hg38


VEP_DEFAULT_PARAMS = {
    "species": "human",
    "variant_class": True,
    "pick": True,
}

VEP_MAX_CHUNK = 200

HTS_COMPRESSED_VALS = ("GZIP", "BGZF")
"""Possible values for the compression field in a :class:`pysam.VariantFile`.
"""

VCF_MASK_COLS = [2, 5, 6, 7]
"""Columns to mask in the VCF file.

Values in these columns are replaced with '.' before querying the
Ensembl API.
"""

API_LONG_RETRY = 60
"""Delay retrying an API call after a 429 response without Retry-After header.

Currently this is the only place in the code where such a delay is defined. In
the future, more settings for the httpx client may be desirable.
"""


@contextlib.contextmanager
def _text_from_bgzf(file_path, mode="r"):
    """Wraps a UTF-8 byte decoder around a BGZF-decompressor.

    This is just a helper to open .vcf.gz files with the same call signature
    as :func:`open` would on a plain .vcf file.
    """
    with pysam.BGZFile(file_path, mode) as bgzf, io.TextIOWrapper(
        bgzf, encoding="utf-8"
    ) as text:
        yield text


def vcf_rows(vcf_path: pathlib.Path) -> typing.Iterator[list[str]]:
    """Yields rows from a VCF file.

    This function opens a VCF file that's either plain text or compressed.
    Normally the :class:`pysam.VariantFile` class would be used to read
    VCF files but in this case there's no need to parse the VCF fully as
    VariantFile does, it's just a matter of reading tab-separated lines
    of text.

    This function could be moved to a more general module in the future,
    but currently it's only needed by the Ensembl API querying code.

    :param vcf_path: Path to the VCF file.
    :return: Iterator of rows from the VCF file.
    """
    file_opener = open
    with pysam.VariantFile(vcf_path, "r") as vf:
        if not vf.is_vcf:
            raise ValueError(f"{vcf_path} does not appear to be a VCF file.")
        if vf.compression in HTS_COMPRESSED_VALS:
            file_opener = _text_from_bgzf
    with file_opener(vcf_path, "r") as vcf:
        reader = csv.reader(vcf, delimiter="\t")
        for row in reader:
            if row and len(row) >= 1 and row[0].startswith("#"):
                continue
            yield row


def vcf_to_vep_query_data(
    vcf_path: pathlib.Path, chunk_size=VEP_MAX_CHUNK
) -> typing.Iterator[list[str]]:
    """Gets lists of data in chunks from the VCF file for querying the Ensembl API.

    This is tailored to provide input to the
    `POST vep/:species/region <https://rest.ensembl.org/documentation/info/vep_region_post>`_
    endpoint.

    :param vcf_path: Path to the VCF file.
    :param chunk_size: Maximum number of variants to include in each chunk.
    """
    chunk = []
    for row in vcf_rows(vcf_path):
        for col in VCF_MASK_COLS:
            row[col] = "."
        chunk.append(" ".join(row[:8]))
        if len(chunk) == chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def query_vep_api(
    client: httpx.Client,
    chunk: list[str],
    assembly: Assembly = Assembly.GRCh37,
    retries=3,
    params: dict | None = VEP_DEFAULT_PARAMS,
    response_extractor: typing.Callable[[dict], dict] | None = None,
):
    """Queries the Ensembl Variant Effect Predictor (VEP) API.

    The API endpoint is `POST vep/:species/region
    <https://rest.ensembl.org/documentation/info/vep_region_post>`_. The API has
    several limits on this endpoint, including a maximum of 200 variants per
    request, and may return a 429 status code if the rate limit is exceeded.  This
    function can retry the request multiple times, using the Retry-After header if
    it's present to delay the next request.

    This function uses the synchronous client from the
    `httpx <https://www.python-httpx.org/>`_ library, which is not as powerful as
    the asynchronous client.  A future version of this function may use the
    asynchronous client to improve performance.

    :param chunk: List of (up to 200) strings from the VCF file.  See
        :func:`get_vcf_query_data` for the format of these strings.
    :param assembly: Assembly to use for the API query, by default GRCh37.
    :param retries: Number of times to retry the API call if it returns a 429 status
        code.
    :param params: Additional parameters to pass to the API.  By default, this
        includes the `species=human` `variant_class=true`, and `pick=true`
        parameters.  `pick=true` is helpful for simplifying the consequences
        of the variants and make it easier to assign a single annotation derived
        from the consequences to the variant.
    :param response_extractor: Optional function to transform data from the API
        response. Without this, there is probably more fields in the response
        with fields without the desired names. Additionally, the default dict
        items returned have additional structure (sublist, subdicts) that is
        preferably flattened by this function.
    """
    subdomain = "grch37.rest" if assembly == Assembly.GRCh37 else "rest"
    url = f"https://{subdomain}.ensembl.org/vep/homo_sapiens/region"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    data = json.dumps({"variants": chunk})
    retried = 0
    while retried <= retries:
        response = client.post(url, headers=headers, data=data, params=params)
        if response.status_code == 200:
            items = response.json()
            if response_extractor:
                # optionally transform data
                items = list(map(response_extractor, items))
            return items
        elif response.status_code == 429:
            retried += 1
            retry_after = response.headers.get("Retry-After")
            if retry_after:
                sleep_time = int(retry_after)
                logger.warn(
                    f"API code 429 with Retry-After. Retrying after {sleep_time} seconds."
                )
                time.sleep(sleep_time)
            else:
                logger.warn(
                    f"API code 429 with no Retry-After. Retrying after {API_LONG_RETRY} seconds."
                )
                time.sleep(sleep_time)
        else:
            response.raise_for_status()
    raise TimeoutError("too many retries")

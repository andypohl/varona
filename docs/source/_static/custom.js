// _static/custom.js
document.addEventListener("DOMContentLoaded", function() {
    var links = document.getElementsByTagName('a');
    for (var i = 0; i < links.length; i++) {
        if (links[i].href.includes("repo.html")) {
            links[i].href = "https://github.com/andypohl/varona";
            links[i].innerHTML = '<i class="fab fa-github"></i> ' + links[i].innerHTML;
        }
    }
});

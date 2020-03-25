var i = 0;
var pos = 0;

function timedCount() {

    postMessage(pos);



    setTimeout("timedCount()", 1);
}

timedCount();
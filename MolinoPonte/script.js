var posx = [];
var posy = [];

var w;
var t;
var totautosx = 1;
var totautodx = 1;

var finito = false;
var finito2 = false;


function RespawnDX() {
    for (i = 0; i < totautodx; i++) {
        var elem = document.getElementById("car" + (i + 3));
        elem.parentNode.removeChild(elem);
    }

    for (i = totautodx - 1; i >= 0; i--) {
        posx[i + 3] = 1900 + i * 60;
        posy[i + 3] = -200;

        var car = document.createElement("img");

        var randdx = Math.floor(Math.random() * 2);

        if (randdx == 0) {
            car.src = "img/3cariso.png";
        }
        else {
            car.src = "img/isocamionDX.png";
            posy[i + 3] = -250;
        }
        car.style = "width: 150px; position: absolute;";
        car.style.transform = "rotate(0deg)";
        car.style.left = posx[i + 3];
        car.style.top = posy[i + 3];


        car.id = "car" + (i + 3);
        document.body.appendChild(car);
    }
}

function RespawnSX() {
    for (i = 0; i < totautosx; i++) {
        var elem = document.getElementById("car" + i);
        elem.parentNode.removeChild(elem);
    }


    for (i = totautosx - 1; i >= 0; i--) {
        posx[i] = -400 + i * 50;
        posy[i] = 1125;

        var car = document.createElement("img");

        var randsx = Math.floor(Math.random() * 2);

        if (randsx == 0) {
            car.src = "img/3cariso2.png";
        }
        else {
            car.src = "img/isocamionSX.png"
            posy[i] = 1125;
        }

        car.style = "width: 150px; position: absolute;";
        car.style.transform = "rotate(0deg)";
        car.style.left = posx[i];
        car.style.top = posy[i];
        car.id = "car" + i;
        document.body.appendChild(car);
    }


}


function newcar() {

    var dotL = document.getElementById("dotL");
    var dotR = document.getElementById("dotR");

    if (typeof (t) == "undefined") {
        t = new Worker("WB_Semaforo.js");
    }
    t.onmessage = function (event) {
        if (event.data == 0) {
            dotL.style.backgroundColor = "green";
            dotR.style.backgroundColor = "red";
        }
        else {
            dotR.style.backgroundColor = "green";
            dotL.style.backgroundColor = "red";
        }

        t.terminate();

    }



    for (i = totautosx - 1; i >= 0; i--) {
        posx[i] = -400 + i * 50;
        posy[i] = 1125;

        var car = document.createElement("img");
        car.src = "img/3cariso2.png";
        car.style = "width: 150px; position: absolute;";
        car.style.transform = "rotate(0deg)";
        car.style.left = posx[i];
        car.style.top = posy[i];
        car.id = "car" + i;

        document.body.appendChild(car);
    }



    for (i = totautodx - 1; i >= 0; i--) {
        posx[i + 3] = 1900 + i * 60;
        posy[i + 3] = -200;

        var car = document.createElement("img");
        car.src = "img/3cariso.png";
        car.style = "width: 150px; position: absolute;";
        car.style.transform = "rotate(0deg)";
        car.style.left = posx[i + 3];
        car.style.top = posy[i + 3];

        car.id = "car" + (i + 3);
        document.body.appendChild(car);
    }

    startWorker();
}



function startWorker() {


    if (typeof (w) == "undefined") {
        w = new Worker("WB_Movement.js");
    }


    posy.fill = 0;


    w.onmessage = function (event) {
        var dotL = document.getElementById("dotL");
        var dotR = document.getElementById("dotR");
        var stop = false;

        var stop3 = false;

        var top0 = document.getElementById("car0").style.top.split("p");
        top0 = Number(top0[0]);

        var left0 = document.getElementById("car0").style.left.split("p");
        left0 = Number(left0[0]);

        var top1 = document.getElementById("car3").style.top.split("p");
        top1 = Number(top1[0]);

        var left1 = document.getElementById("car3").style.left.split("p");
        left1 = Number(left1[0]);

        if (dotL.style.backgroundColor == "green")

            if (top0 <= -400 && !stop) {
                stop = true;

                finito = true;
            }

        if (dotR.style.backgroundColor == "green")
            if (top1 >= 1000 && !stop3) {
                stop3 = true;
                finito2 = true;
            }



        if (dotL.style.backgroundColor == "green") {
            if (!stop) {
                for (var k = 0; k < totautosx; k++) {
                    posy[k] -= 0.29 * 2;
                    posx[k] += 0.5 * 2;
                    var car = document.getElementById("car" + k);
                    car.style.top = posy[k];
                    car.style.left = posx[k];
                }
            }

        }

        if (dotR.style.backgroundColor == "green") {
            if (!stop3) {
                for (var k = 0; k < totautodx; k++) {
                    posy[k + 3] += 0.29 * 2;
                    posx[k + 3] -= 0.5 * 2;
                    var car = document.getElementById("car" + (k + 3));
                    car.style.top = posy[k + 3];
                    car.style.left = posx[k + 3];

                }
            }

        }



        if (finito) {
            dotL.style.backgroundColor = "red";
            dotR.style.backgroundColor = "green";
            RespawnSX();
            finito = false;
        }
        if (finito2) {
            dotL.style.backgroundColor = "green";
            dotR.style.backgroundColor = "red";
            RespawnDX();
            finito2 = false;
        }
    }




}








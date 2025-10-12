const canvas = document.getElementById('puyogrid');
const ctx = canvas.getContext('2d');
const width = canvas.width = 800;
const height = canvas.height = 600;
const container = { left: 100, right: 700, top: 0, bottom: height };
const r = 15; // ball radius
let balls = [];
let N = 6;
const colors = ["transparent", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#8800FF", "#88FF00", "#FF0088", "#FF8800", "#00FF88", "#0088FF", "#888800", "#880088", "#008888", "#FF8888", "#8888FF", "#88FF88"];

// Matter.js setup
const engine = Matter.Engine.create();
engine.world.gravity.y = 2; // Stronger gravity
const world = engine.world;

// Create walls aligned with container
const wallOptions = { isStatic: true, render: { visible: false } };
const wallThickness = 10;
const leftWall = Matter.Bodies.rectangle(container.left + wallThickness / 2, height / 2, wallThickness, height, wallOptions);
const rightWall = Matter.Bodies.rectangle(container.right - wallThickness / 2, height / 2, wallThickness, height, wallOptions);
const bottomWall = Matter.Bodies.rectangle(width / 2, container.bottom - wallThickness / 2, width, wallThickness, wallOptions);
Matter.World.add(world, [leftWall, rightWall, bottomWall]);

const controls = document.getElementById('controls');
controls.innerHTML = '';
const label = document.createElement('label');
label.textContent = 'N:';
label.style.fontFamily = 'PuyoFont';
label.style.padding = '4px';
label.style.marginRight = '10px';
controls.appendChild(label);
const nInput = document.createElement('input');
nInput.type = 'number';
nInput.value = N;
nInput.min = 1;
nInput.max = colors.length - 1;
nInput.style.width = '3.5rem';
nInput.style.fontSize = '2rem';
nInput.style.fontFamily = 'PuyoFont';
nInput.style.backgroundColor = '#222222';
nInput.style.color = '#fac9d5';
nInput.style.borderRight = '2px solid #444444';
controls.appendChild(nInput);
const resetBtn = document.createElement('button');
resetBtn.textContent = 'Reset';
resetBtn.style.color = '#fac9d5';
resetBtn.style.marginLeft = '10px';
resetBtn.addEventListener('click', () => {
    N = parseInt(nInput.value);
    if (isNaN(N) || N < 1 || N > colors.length - 1) {
        alert('Invalid N');
        return;
    }
    balls = [];
});
controls.appendChild(resetBtn);

function draw() {
    ctx.clearRect(0, 0, width, height);
    // draw container walls (left, right, bottom)
    ctx.strokeStyle = 'white';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(container.left, container.top);
    ctx.lineTo(container.left, container.bottom);
    ctx.moveTo(container.right, container.top);
    ctx.lineTo(container.right, container.bottom);
    ctx.moveTo(container.left, container.bottom);
    ctx.lineTo(container.right, container.bottom);
    ctx.stroke();
    // draw balls
    balls.forEach(ball => {
        if (ball.state === 'removed') return;
        const body = ball.body;
        ctx.beginPath();
        ctx.arc(body.position.x, body.position.y, r, 0, 2 * Math.PI);
        let fillColor = colors[ball.color];
        if (ball.state === 'eliminating') {
            const progress = Math.min((Date.now() - ball.elimTime) / 500, 1);
            const origColor = hexToRgb(colors[ball.color]);
            const white = { r: 255, g: 255, b: 255 };
            const rInterp = Math.round(origColor.r + (white.r - origColor.r) * progress);
            const gInterp = Math.round(origColor.g + (white.g - origColor.g) * progress);
            const bInterp = Math.round(origColor.b + (white.b - origColor.b) * progress);
            fillColor = `rgb(${rInterp}, ${gInterp}, ${bInterp})`;
        }
        ctx.fillStyle = fillColor;
        ctx.fill();
    });
}

let lastTime = 0;
function update(currentTime) {
    const deltaTime = currentTime - lastTime;
    lastTime = currentTime;
    Matter.Engine.update(engine, deltaTime);

    // check for elimination
    for (let i = 0; i < balls.length; i++) {
        for (let j = i + 1; j < balls.length; j++) {
            const b1 = balls[i];
            const b2 = balls[j];
            if (b1.state === 'normal' && b2.state === 'normal' && b1.color === b2.color) {
                const dx = b2.body.position.x - b1.body.position.x;
                const dy = b2.body.position.y - b1.body.position.y;
                const dist = Math.sqrt(dx * dx + dy * dy);
                if (dist <= 2 * r + 1) { // touching
                    b1.state = 'eliminating';
                    b2.state = 'eliminating';
                    b1.elimTime = Date.now();
                    b2.elimTime = Date.now();
                }
            }
        }
    }

    // update eliminating
    const now = Date.now();
    balls.forEach(ball => {
        if (ball.state === 'eliminating' && now - ball.elimTime >= 500) {
            ball.state = 'removed';
            Matter.World.remove(world, ball.body); // Remove from physics
        }
    });

    // remove eliminated balls
    balls = balls.filter(ball => ball.state !== 'removed');
}

function loop(currentTime) {
    update(currentTime);
    draw();
    requestAnimationFrame(loop);
}

document.addEventListener('keydown', (event) => {
    if (event.code === 'Space') {
        const x = container.left + r + Math.random() * (container.right - container.left - 2 * r);
        const y = container.top - r;
        const color = Math.floor(Math.random() * N) + 1;
        const body = Matter.Bodies.circle(x, y, r, {
            restitution: 0.1,
            frictionAir: 0.02,
            density: 0.001
        });
        Matter.World.add(world, body);
        balls.push({ body, color, state: 'normal' });
    }
    if (event.code === 'ArrowDown') {
        N = Math.max(1, N - 1);
        nInput.value = N;
    }
    if (event.code === 'ArrowUp') {
        N = Math.min(colors.length - 1, N + 1);  // Max at 15
        nInput.value = N;
    }
    if (event.code === 'KeyR') {
        balls = [];
    }
});

requestAnimationFrame(loop);

function hexToRgb(hex) {
    const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? {
        r: parseInt(result[1], 16),
        g: parseInt(result[2], 16),
        b: parseInt(result[3], 16)
    } : null;
}
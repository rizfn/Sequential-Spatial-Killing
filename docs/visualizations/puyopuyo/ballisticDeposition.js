// Create a 20x20 grid
const L = 20;
let N = 6;
const grid = Array.from({ length: L }, () => Array(L).fill(0));
const colors = ["transparent", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#8800FF", "#FF8800", "#00FF88", "#888800", "#8800FF", "#0088FF", "#FF0088", "#8888FF", "#FF8888"];

// Dynamically calculate cell size based on window dimensions
const cellSize = Math.min(window.innerWidth / L, window.innerHeight / L);

const svg = d3.select("svg#puyogrid");

const gridGroup = svg.append("g")
    .attr("transform", `translate(${(window.innerWidth - L * cellSize) / 2}, ${(window.innerHeight - L * cellSize) / 2})`);

const controls = d3.select("div#controls");

controls.append("label")
    .text("N:")
    .style("font-family", "PuyoFont")
    .style("padding", "4px")
    .style("margin-right", "10px");

const nInput = controls.append("input")
    .attr("type", "number")
    .attr("value", N)
    .attr("min", 1)
    .attr("max", colors.length - 1)
    .style("width", "3.5rem")
    .style("font-size", "2rem")
    .style("font-family", "PuyoFont")
    .style("background-color", "#222222")
    .style("color", "#fac9d5")
    .style("border-right", "2px solid #444444");

controls.append("button")
    .text("Reset")
    .style("color", "#fac9d5")
    .style("margin-left", "10px")
    .on("click", () => {
        N = parseInt(nInput.property("value"), 10);
        if (isNaN(N) || N < 1 || N > colors.length - 1) {
            alert("Please enter a valid number for N (between 1 and " + (colors.length - 1) + ").");
            return;
        }
        grid.forEach(row => row.fill(0));
        currentPuyo = null;
        if (animationFrame) {
            cancelAnimationFrame(animationFrame);
            animationFrame = null;
        }
        drawGrid();
    });


function drawGrid() {
    gridGroup.selectAll("rect")
        .data(grid.flat())
        .join("rect")
        .attr("x", (_, i) => (i % L) * cellSize)
        .attr("y", (_, i) => Math.floor(i / L) * cellSize)
        .attr("width", cellSize)
        .attr("height", cellSize)
        .attr("fill", d => colors[d]);
}

// Find connected components of the same species
function connectedComponents(grid) {
    const labels = Array.from({ length: L }, () => Array(L).fill(0));
    let label = 0;

    function dfs(row, col, currentLabel, species) {
        if (
            row < 0 || row >= L || col < 0 || col >= L ||
            grid[row][col] !== species || labels[row][col] !== 0
        ) {
            return;
        }
        labels[row][col] = currentLabel;
        dfs(row - 1, col, currentLabel, species);
        dfs(row + 1, col, currentLabel, species);
        dfs(row, col - 1, currentLabel, species);
        dfs(row, col + 1, currentLabel, species);
    }

    for (let row = 0; row < L; row++) {
        for (let col = 0; col < L; col++) {
            if (grid[row][col] !== 0 && labels[row][col] === 0) {
                label++;
                dfs(row, col, label, grid[row][col]);
            }
        }
    }

    return { labels, N: label };
}

function placePuyo(column, species) {
    if (grid[0][column] !== 0) return false; // Column is full
    grid[0][column] = species;
    return true;
}

function removePuyo() {
    const { labels, N } = connectedComponents(grid);
    for (let label = 1; label <= N; label++) {
        const cells = [];
        for (let row = 0; row < L; row++) {
            for (let col = 0; col < L; col++) {
                if (labels[row][col] === label) {
                    cells.push([row, col]);
                }
            }
        }
        if (cells.length > 1) {
            const species = grid[cells[0][0]][cells[0][1]];
            if (cells.some(([r, c]) => grid[r][c] === species)) {
                // Animate the removal
                cells.forEach(([row, col]) => {
                    const rect = gridGroup.select(`rect[x="${col * cellSize}"][y="${row * cellSize}"]`);
                    rect.transition()
                        .duration(300)
                        .attr("width", cellSize * 1.5)
                        .attr("height", cellSize * 1.5)
                        .attr("x", col * cellSize - cellSize * 0.25)
                        .attr("y", row * cellSize - cellSize * 0.25)
                        .attr("fill", "#FFFFFF")
                        .on("end", () => {
                            grid[row][col] = 0; // Set the cell to empty after animation
                            drawGrid(); // Redraw the grid
                        });
                });
            }
        }
    }
}

function fall() {
    let changed = false;
    for (let row = L - 2; row >= 0; row--) { // Start from the second-to-last row
        for (let col = 0; col < L; col++) {
            if (grid[row][col] > 0) { // If there's a cell at this position
                const species = grid[row][col];

                // Check if the cell can fall
                if (
                    grid[row + 1][col] === 0 && // Below is empty
                    (col === 0 || grid[row][col - 1] === 0) && // No neighbor to the bottom-left
                    (col === L - 1 || grid[row][col + 1] === 0) && // No neighbor to the bottom-right
                    (col === 0 || grid[row][col - 1] === 0) && // No neighbor to the left
                    (col === L - 1 || grid[row][col + 1] === 0) // No neighbor to the right
                ) {
                    // Move the cell down
                    grid[row + 1][col] = species;
                    grid[row][col] = 0;
                    changed = true;
                }
            }
        }
    }
    return changed;
}

let currentPuyo = null;
let animationFrame = null;

function step() {
    if (currentPuyo) {
        const { column, species, row } = currentPuyo;

        // Check if the cell can fall further
        if (
            row < L - 1 &&
            grid[row + 1][column] === 0 && // Check below
            (column === 0 || grid[row][column - 1] === 0) && // No neighbor to the bottom-left
            (column === L - 1 || grid[row][column + 1] === 0) && // No neighbor to the bottom-right
            (column === 0 || grid[row][column - 1] === 0) && // No neighbor to the left
            (column === L - 1 || grid[row][column + 1] === 0) // No neighbor to the right
        ) {
            grid[row][column] = 0;
            grid[row + 1][column] = species;
            currentPuyo.row++;
        } else {
            // Stop falling and reset currentPuyo
            currentPuyo = null;
            removePuyo();
        }
    } else {
        // Apply the fall function to all cells
        if (!fall()) {
            cancelAnimationFrame(animationFrame);
            return;
        }
    }

    drawGrid();
    animationFrame = requestAnimationFrame(step);
}

document.addEventListener("keydown", (event) => {
    if (event.code === "Space" && !currentPuyo) {
        const column = Math.floor(Math.random() * L);
        const species = Math.floor(Math.random() * N) + 1;
        if (placePuyo(column, species)) {
            currentPuyo = { column, species, row: 0 };
            animationFrame = requestAnimationFrame(step);
        }
    }
    // on down arrow, decrease N, on up arrow, increase N
    if (event.code === "ArrowDown") {
        N = Math.max(1, N - 1);
        nInput.property("value", N);
    }
    if (event.code === "ArrowUp") {
        N = Math.min(colors.length - 1, N + 1);
        nInput.property("value", N);
    }
    // on R, reset sim
    if (event.code === "KeyR") {
        grid.forEach(row => row.fill(0));
        currentPuyo = null;
        if (animationFrame) {
            cancelAnimationFrame(animationFrame);
            animationFrame = null;
        }
        drawGrid();
    }
});


// Initial draw
drawGrid();
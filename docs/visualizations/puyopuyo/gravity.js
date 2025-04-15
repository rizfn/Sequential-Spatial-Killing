// Create a 20x20 grid
const L = 20;
const grid = Array.from({ length: L }, () => Array(L).fill(0));
const colors = ["transparent", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF"]; // Colors for species

// Dynamically calculate cell size based on window dimensions
const cellSize = Math.min(window.innerWidth / L, window.innerHeight / L);

// Initialize D3.js SVG
const svg = d3.select("body")
    .append("svg")
    .attr("width", window.innerWidth)
    .attr("height", window.innerHeight)
    .style("display", "flex")
    .style("align-items", "center")
    .style("justify-content", "center")
    .style("background-image", "linear-gradient(rgba(0, 0, 0, 0.8), rgba(0, 0, 0, 0.8)), url('img/Witch_PPT2.png')")
    .style("background-size", "contain")
    .style("background-repeat", "no-repeat")
    .style("background-position", "center");

    const gridGroup = svg.append("g")
    .attr("transform", `translate(${(window.innerWidth - L * cellSize) / 2}, ${(window.innerHeight - L * cellSize) / 2})`);

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
    for (let row = L - 2; row >= 0; row--) {
        for (let col = 0; col < L; col++) {
            if (grid[row][col] > 0 && grid[row + 1][col] === 0) {
                grid[row + 1][col] = grid[row][col];
                grid[row][col] = 0;
                changed = true;
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
        if (row < L - 1 && grid[row + 1][column] === 0) {
            grid[row][column] = 0;
            grid[row + 1][column] = species;
            currentPuyo.row++;
        } else {
            currentPuyo = null;
            removePuyo();
        }
    } else {
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
        const species = Math.floor(Math.random() * (colors.length - 1)) + 1;
        if (placePuyo(column, species)) {
            currentPuyo = { column, species, row: 0 };
            animationFrame = requestAnimationFrame(step);
        }
    }
});

// Initial draw
drawGrid();
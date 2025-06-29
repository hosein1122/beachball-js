import { beach, plotDC } from "../src/beachball.js";
import { MomentTensor, NodalPlane } from "../src/classes.js";
import { writeFileSync } from "fs";
import { promises as fs } from "fs";
import { createCanvas } from "canvas";

// const Mrr = -58140500000000000;
// const Mtt = 74811200000000000;
// const Mpp = -16670600000000000;
// const Mrt = -7755400000000000;
// const Mrp = 16777100000000000;
// const Mtp = 16055700000000000;

const strike = 240;
const dip = 46;
const rake = -119;

const np1 = new NodalPlane(strike, dip, rake);
const dc = plotDC(np1);

// console.log(dc);

// console.log('first 5 verts tension', dc.tension.vertices.slice(0,5));


writeFileSync("dc_js.json", JSON.stringify(dc));
console.log("✅ dc_js.json written");


/* ───── Helper: read JSON file and return [colors, patches] ───── */
async function loadBeachball(path) {
    const raw = await fs.readFile(path, "utf8");
    const [colors, patches] = JSON.parse(raw);
    return { colors, patches };
}

/* ───── Helper: draw a single patch (closed polygon) ───── */
function drawPatch(ctx, patch, color, offsetX, offsetY) {
    const { vertices } = patch;

    ctx.beginPath();
    /*  Coordinate-system note
    ──────────────────────
    • Python / Matplotlib:
        – Origin (0, 0) is at the BOTTOM-LEFT of the axes.
        – The Y-axis points UP (positive values go upward).

    • HTML5 Canvas (browser or Node “canvas” module):
        – Origin (0, 0) is at the TOP-LEFT of the canvas element.
        – The Y-axis points DOWN (positive values go downward).

    Result: a figure generated with Matplotlib appears vertically flipped
    when drawn directly on a Canvas.  To match the Python view you must
    either (a) invert Y for every vertex (e.g.  ctx.lineTo(x, offsetY - y))
    or (b) apply one global transform:

          ctx.translate(0, canvas.height);   // move origin to bottom-left
          ctx.scale(1, -1);                  // invert Y-axis

    After that transform you can use the same (x, y) coordinates that work
    in Matplotlib, but remember to restore or re-invert the scale before
    drawing text or UI elements that should not be flipped.
    */
    // Y is inverted:  offsetY - y
    ctx.moveTo(vertices[0][0] + offsetX, offsetY - vertices[0][1]);   // ← changed
    for (let i = 1; i < vertices.length; i++) {
        ctx.lineTo(vertices[i][0] + offsetX, offsetY - vertices[i][1]); // ← changed
    }
    ctx.closePath();

    ctx.fillStyle = color === "b" ? "#0000ff" : "#ffffff";
    ctx.fill();
    ctx.lineWidth = 0.5;
    ctx.strokeStyle = "#000000";
    ctx.stroke();
}


/* ───── Main routine ───── */
(async () => {
    // 1) load both JSON files
    const py = await loadBeachball("dc_py.json");
    const js = await loadBeachball("dc_js.json");

    // 2) canvas geometry
    const SIZE = 220;           // width/height for each beach-ball
    const PAD = 20;            // spacing between them
    const W = SIZE * 2 + PAD * 3;
    const H = SIZE + PAD * 2;

    const canvas = createCanvas(W, H);
    const ctx = canvas.getContext("2d");

    ctx.fillStyle = "#eeeeee";  // light background
    ctx.fillRect(0, 0, W, H);

    // 3) draw Python beach-ball (left)
    for (let i = 0; i < py.patches.length; i++) {
        drawPatch(
            ctx,
            py.patches[i],
            py.colors[i],
            PAD + SIZE / 2,          // center X
            PAD + SIZE / 2           // center Y
        );
    }

    // 4) draw JavaScript beach-ball (right)
    for (let i = 0; i < js.patches.length; i++) {
        drawPatch(
            ctx,
            js.patches[i],
            js.colors[i],
            PAD * 2 + SIZE + SIZE / 2,
            PAD + SIZE / 2
        );
    }

    // 5) annotate titles
    ctx.fillStyle = "#000000";
    ctx.font = "16px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Python", PAD + SIZE / 2, PAD / 2 + 16);
    ctx.fillText("JavaScript", PAD * 2 + SIZE + SIZE / 2, PAD / 2 + 16);

    // 6) save as PNG
    const out = canvas.toBuffer("image/png");
    await fs.writeFile("compare.png", out);
    console.log("✅  compare.png written.");
})();


// const mt = new MomentTensor(Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, 0);
// const col1 = beach(mt, { radius: 100 });

// console.log(col1);

// writeFileSync("col1.json", JSON.stringify(col1));
// console.log("✅   col1.json written");
// file: test3.js
import { beach, plotDC, plotMT, mt2axes } from "../src/beachball.js";
import { MomentTensor, NodalPlane }   from "../src/classes.js";
import fs                             from "fs";
import { createCanvas } from "canvas";

const Mrr = -58140500000000000;
const Mtt =  74811200000000000;
const Mpp = -16670600000000000;
const Mrt =  -7755400000000000;
const Mrp =  16777100000000000;
const Mtp =  16055700000000000;

// 1) بسازیم MomentTensor
const mt = new MomentTensor(Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, 0);

// 2) بگیریم خروجی کلی beach
const bch = beach(mt);
// console.log(bch);
console.log(bch.patches[1].vertices);


// 3) ذخیره‌سازی JSON
fs.writeFileSync(
  "beach_test3.json",
  JSON.stringify(bch, null, 2),
  "utf-8"
);
console.log("✅ beach_test3.json saved");


// file: draw_test3.js
const WIDTH  = 400;   // ابعاد بوم
const HEIGHT = 400;
const CENTER = [WIDTH/2, HEIGHT/2];  // مرکز ترسیم

// 1) بارگذاری JSON
const data = JSON.parse(fs.readFileSync("beach_test3.json", "utf-8"));
const { patches, fillColors, edgecolor, linewidth, alpha } = data;

// 2) آماده‌سازی Canvas
const canvas = createCanvas(WIDTH, HEIGHT);
const ctx    = canvas.getContext("2d");

// پس‌زمینه
ctx.fillStyle = "white";
ctx.fillRect(0, 0, WIDTH, HEIGHT);

// 3) تابع عمومی برای رسم پچ‌ها
function drawPatches(ctx) {
  ctx.globalAlpha = alpha;
  ctx.lineWidth   = linewidth;
  ctx.strokeStyle = edgecolor;

  for (let i = 0; i < patches.length; i++) {
    const p   = patches[i];
    const fc  = fillColors[i];
    if (fc) ctx.fillStyle = fc;

    ctx.beginPath();
    if (p.type === "ellipse") {
      // ⚠ Note: canvas Y grows downward; matplotlib Y grows upward
      ctx.ellipse(
        CENTER[0] + p.xy[0],
        CENTER[1] - p.xy[1],
        p.width  / 2,
        p.height / 2,
        0, 0, 2*Math.PI
      );
      if (fc) ctx.fill();
      ctx.stroke();
    } else if (p.vertices) {
      const verts = p.vertices;
      // شروع مسیر
      let [x0,y0] = verts[0];
      ctx.moveTo(
        CENTER[0] + x0 * 1, 
        CENTER[1] - y0 * 1
      );
      for (let j = 1; j < verts.length; j++) {
        let [x, y] = verts[j];
        ctx.lineTo(
          CENTER[0] + x,
          CENTER[1] - y
        );
      }
      ctx.closePath();
      if (fc) ctx.fill();
      ctx.stroke();
    }
  }
}

// 4) رسم و ذخیره PNG
drawPatches(ctx);

const out = fs.createWriteStream("beach_test3.png");
canvas.createPNGStream().pipe(out);
out.on("finish", () => console.log("✅ beach_test3.png written"));


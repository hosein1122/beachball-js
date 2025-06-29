import { beach, plotDC, plotMT, mt2axes } from "../src/beachball.js";
import { MomentTensor, NodalPlane } from "../src/classes.js";
import fs from 'fs';
import { createCanvas } from "canvas";

const Mrr = -58140500000000000;
const Mtt = 74811200000000000;
const Mpp = -16670600000000000;
const Mrt = -7755400000000000;
const Mrp = 16777100000000000;
const Mtp = 16055700000000000;

const strike = 240;
const dip = 46;
const rake = -119;

// const np1 = new NodalPlane(strike, dip, rake);
// const dc = plotDC(np1);

// console.log(dc);

// console.log('first 5 verts tension', dc.tension.vertices.slice(0,5));



const mt = new MomentTensor(Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, 0);
const [T, N, P] = mt2axes(mt.normalized);

const pmt = plotMT(T, N, P);
// console.log(pmt);

// ذخیره در فایل JSON
fs.writeFileSync("pmt.json", JSON.stringify(pmt, null, 2));
console.log("✅ Saved to pmt.json");





// const { createCanvas } = require('canvas');
// const fs = require('fs');

const canvas = createCanvas(500, 300); // عرض، ارتفاع
const ctx = canvas.getContext('2d');

// Helper: ترسیم بیچ‌بال
function drawBeachball(data, centerX, centerY) {
  const [colors, patches] = data;

  for (let i = 0; i < patches.length; i++) {
    const patch = patches[i];
    const color = colors[i];
    ctx.beginPath();
    ctx.fillStyle = color === 'b' ? 'black' : 'white';
    ctx.strokeStyle = 'black';

    if (patch.type === 'ellipse') {
      ctx.ellipse(
        centerX,
        centerY,
        patch.width / 2,
        patch.height / 2,
        0, 0, 2 * Math.PI
      );
      ctx.fill();
      ctx.stroke();
    } else {
      const verts = patch.vertices;
      if (verts.length === 0) continue;
      const [x0, y0] = verts[0];
      ctx.moveTo(centerX + x0 - patch.center[0], centerY - (y0 - patch.center[1]));

      for (let j = 1; j < verts.length; j++) {
        const [x, y] = verts[j];
        ctx.lineTo(centerX + x - patch.center[0], centerY - (y - patch.center[1]));
      }
      ctx.closePath();
      ctx.fill();
      ctx.stroke();
    }
  }
}

// Load both JSON files and draw them
const dataPy = JSON.parse(fs.readFileSync('pmt_py.json', 'utf8'));
const dataJs = JSON.parse(fs.readFileSync('pmt.json', 'utf8'));

// رسم هر دو بیچ‌بال
drawBeachball(dataPy, 125, 150);  // سمت چپ
drawBeachball(dataJs, 375, 150);  // سمت راست

// ذخیره تصویر
const out = fs.createWriteStream('beachball_compare.png');
const stream = canvas.createPNGStream();
stream.pipe(out);
out.on('finish', () => console.log('Saved as beachball_compare.png'));

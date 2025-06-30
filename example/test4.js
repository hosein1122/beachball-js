// render.js
import { createCanvas } from 'canvas';
import fs from 'fs';
import { beach } from '../src/beachball.js';
import { MomentTensor } from '../src/classes.js';

// ─── 1) تابع رسم beachball ───────────────────────────────────────────────
/**
 * رسم یک patch (ellipse یا polygon) روی کانواس
 * بر مبنای ساختار خروجیِ beach().patches[i]
 *
 * @param {CanvasRenderingContext2D} ctx
 * @param {object}                patch    – { type:'ellipse', xy, width, height } | { vertices: [...], ... }
 * @param {'b'|'w'|null}          fillFlag – از data.fillColors[i]
 * @param {number}                offsetX  – مقدار جا‌به‌جایی X (مثلاً data.xy[0])
 * @param {number}                offsetY  – مقدار جا‌به‌جایی Y (مثلاً data.xy[1])
 */
function drawPatch(ctx, patch, fillFlag, offsetX, offsetY) {
    ctx.beginPath();

    if (patch.type === 'ellipse') {
        // بیضی
        const [cx, cy] = patch.xy;
        ctx.ellipse(
            cx + offsetX,
            cy + offsetY,
            patch.width,
            patch.height,
            0, 0, 2 * Math.PI
        );

    } else if (patch.vertices) {
        // چندضلعی
        const verts = patch.vertices;
        ctx.moveTo(verts[0][0] + offsetX, verts[0][1] + offsetY);
        for (let i = 1; i < verts.length; i++) {
            ctx.lineTo(verts[i][0] + offsetX, verts[i][1] + offsetY);
        }
        ctx.closePath();
    }

    // پرکردن (b=آبی، w=سفید)
    if (fillFlag === 'b' || fillFlag === 'w') {
        ctx.fillStyle = (fillFlag === 'b' ? '#0000ff' : '#ffffff');
        ctx.fill();
    }

    // رسم لبه با ضخامت و رنگ ثابت
    ctx.lineWidth = 0.5;
    ctx.strokeStyle = '#000000';
    ctx.stroke();
}

/**
 * رسم کامل beachball از آبجکت data = beach(...)
 *
 * @param {CanvasRenderingContext2D} ctx
 * @param {object} data – { patches, fillColors, xy, size, … }
 */
function drawBeach(ctx, data) {
    const { patches, fillColors, xy, size } = data;

    // ۱) یک بار invert محور Y و انتقال origin به پایین-چپ
    ctx.save();
    ctx.translate(0, size);
    ctx.scale(1, -1);

    // ۲) برای هر پچ فراخوانی drawPatch با offset از data.xy
    const [dx, dy] = xy;
    patches.forEach((patch, i) => {
        drawPatch(ctx, patch, fillColors[i], dx, dy);
    });

    ctx.restore();
}



// ─── 2) بخش اجرایی: تولید داده + ذخیره به PNG ───────────────────────────
(async () => {
    // مقادیر Mxx را مطابق درخواست شما اینجا قرار دهید:
    const Mrr = -58140500000000000;
    const Mtt = 74811200000000000;
    const Mpp = -16670600000000000;
    const Mrt = -7755400000000000;
    const Mrp = 16777100000000000;
    const Mtp = 16055700000000000;

    const mt = new MomentTensor(Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, 0);
    const data = beach(mt);

    console.log(data);
    // console.log(data.patches[1].vertices);

    // اندازهٔ دلخواه کانواس (مثلاً 400px)
    const canvas = renderBeach(data, 400);

    // ذخیرهٔ تصویر
    fs.writeFileSync('beachball.png', canvas.toBuffer('image/png'));
    console.log('✅ Beachball saved as beachball.png');
})();

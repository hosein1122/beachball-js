/**
 * Converts patch object { vertices, codes } to a Path2D instance.
 * @param {{vertices: number[][], codes: string[]}} patch
 * @returns {Path2D}
 */
export function patchToPath2D(patch) {
    const path = new Path2D();
    patch.vertices.forEach(([vx, vy], idx) => {
        const code = patch.codes[idx];
        if (code === 'M') path.moveTo(vx, vy);
        else if (code === 'L') path.lineTo(vx, vy);
        else if (code === 'Z') path.closePath();
    });
    return path;
}



//و هنگام رسم در مرورگر:

// import { xy2patch } from './beachball.js';
// import { patchToPath2D } from './canvasUtils.js';

// const patch = xy2patch(x, y, [sx, sy], [cx, cy]);
// const path2d = patchToPath2D(patch);

// ctx.fillStyle = 'blue';
// ctx.fill(path2d);
// ctx.stroke(path2d);

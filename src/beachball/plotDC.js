
// plotDC.js :

import { auxPlane, xy2patch } from "../math/geometry.js";
import { D2R, pol2cart } from "../math/trigHelpers.js";


/**
 * Build the geometric description of a double-couple focal-mechanism
 * “beach-ball” suitable for later rendering on any 2-D backend
 * (Canvas, SVG, WebGL …).
 *
 * No drawing is performed here; the function returns *data only*.
 *
 * @param {Object}        np1              – First nodal-plane object
 * @param {number}        np1.strike       – Strike  (degrees, 0–360)
 * @param {number}        np1.dip          – Dip     (degrees, 0–90)
 * @param {number}        np1.rake         – Rake    (degrees, −180–180 or 0–360)
 * @param {number} [size = 200]            – Output square size in pixels
 * @param {[number,number]} [xy = [0,0]]   – Centre of the plot in user units
 * @param {number|[number,number]} [width = 200]
 *        If a single number, circle radius in pixels.
 *        If a pair [xw, yw], ellipse half-widths along *x* and *y*.
 *
 * @returns {[string[], Object[]]}
 *   • **colours** — two-element array `['b','w']` (pressure colour first).  
 *   • **patches** — array of objects produced by `xy2patch()`, each with:
 *     ─ `vertices` : `[x,y][]` closed polygon (screen units)  
 *     ─ `res`      : `[rx,ry]` per-axis resolution scaling (unitless)  
 *     ─ `center`   : `[cx,cy]` translation applied when creating vertices
 */
export function plotDC(np1, size = 200, xy = [0, 0], width = 200) {
    // ─── 1. Check whether 'width' indicates a Circle (one value) or an Ellipse (two values) ───
    try {
        // If 'width' is a single primitive (e.g. number), it has no 'length' property
        if (width.length === undefined) {
            throw new TypeError();          // Mirrors Python’s TypeError in the assert
        }

        // If iterable, it must contain exactly two elements
        if (width.length !== 2) {
            throw new Error('width must contain exactly two elements'); // Like AssertionError
        }
    } catch (err) {
        // Handle only the TypeError that signals a single-value width
        if (err instanceof TypeError) {
            width = [width, width];         // Same width in both directions → circle
        } else {
            throw err;                      // Re-throw any other errors
        }
    }

    // ─── 2. Normalise rake angle to the 0–180° interval ───
    let s_1 = np1.strike;       // Strike
    let d_1 = np1.dip;          // Dip
    let r_1 = np1.rake;         // Rake

    let m = 0;                  // Flag: 1 if rake was adjusted

    // If rake exceeds 180°, wrap it down into the 0–180° range
    if (r_1 > 180) {
        r_1 -= 180;
        m = 1;                    // Original rake was out of range
    }

    // If rake is negative, shift it up into the 0–180° range
    if (r_1 < 0) {
        r_1 += 180;
        m = 1;                    // Original rake was out of range
    }


    // Get azimuth and dip of the second plane
    let [s_2, d_2, _r_2] = auxPlane(s_1, d_1, r_1); 

    // ─── 3. Prepare half-size and clamp dip angles below 90° ───
    const d = size / 2;                 // Half of the canvas size (radius in pixels)

    // Dip angles of exactly 90° cause singularities; cap them at 89.9999°
    if (d_1 >= 90) d_1 = 89.9999;
    if (d_2 >= 90) d_2 = 89.9999;


    // ─── 4. Create a phi array (0 … π, step 0.01) ───
    // Note: Math.PI is not an exact multiple of 0.01, so the last value is < π
    const phi = [];
    for (let t = 0; t < Math.PI; t += 0.01) {
        phi.push(t);
    }

    // ─── 5. Compute l1 and l2 for the two nodal planes ───
    const ninetySquared = 90 * 90;
    const k1 = 90 - d_1;               // Constant term for first plane
    const k2 = 90 - d_2;               // Constant term for second plane

    const l1 = phi.map(t => {
        // numerator = (90 − d1)²
        // denominator = sin²φ + cos²φ · (90 − d1)² / 90²
        const numerator = k1 * k1;
        const denominator =
            Math.sin(t) ** 2 +
            (Math.cos(t) ** 2) * (k1 * k1) / ninetySquared;
        return Math.sqrt(numerator / denominator);
    });

    const l2 = phi.map(t => {
        const numerator = k2 * k2;
        const denominator =
            Math.sin(t) ** 2 +
            (Math.cos(t) ** 2) * (k2 * k2) / ninetySquared;
        return Math.sqrt(numerator / denominator);
    });

    // ─── 6. Build polygon paths for tension (white) and pressure (black) lobes ───
    const collect = [];

    // Iterate first over the tension side, then over the pressure side
    for (const m_ of [(m + 1) % 2, m]) {
        let inc = 1;

        // Segment on first nodal great-circle
        const [x_1, y_1] = pol2cart(
            phi.map(t => t + s_1 * D2R),   // φ shifted by strike of plane-1
            l1                              // corresponding radii
        );

        let th1, th2, xs_1, ys_1, x_2, y_2;

        if (m_ === 1) {
            // ----- Tension lobe -----
            let lo = s_1 - 180;
            let hi = s_2;
            if (lo > hi) inc = -1;

            th1 = range(s_1 - 180, s_2, inc);           // mimic np.arange
            [xs_1, ys_1] = pol2cart(
                th1.map(t => t * D2R),
                Array(th1.length).fill(90)
            );

            [x_2, y_2] = pol2cart(
                phi.map(t => t + s_2 * D2R),
                l2
            );

            th2 = range(s_2 + 180, s_1, -inc);
        } else {
            // ----- Pressure lobe -----
            let hi = s_1 - 180;
            let lo = s_2 - 180;
            if (lo > hi) inc = -1;

            th1 = range(hi, lo, -inc);
            [xs_1, ys_1] = pol2cart(
                th1.map(t => t * D2R),
                Array(th1.length).fill(90)
            );

            [x_2, y_2] = pol2cart(
                phi.map(t => t + s_2 * D2R),
                l2
            );
            x_2.reverse();                     // reverse to preserve winding
            y_2.reverse();

            th2 = range(s_2, s_1, inc);
        }

        // Closing arc on the auxiliary great-circle
        const [xs_2, ys_2] = pol2cart(
            th2.map(t => t * D2R),
            Array(th2.length).fill(90)
        );

        // Concatenate all pieces into one polygon
        const x = [...x_1, ...xs_1, ...x_2, ...xs_2];
        const y = [...y_1, ...ys_1, ...y_2, ...ys_2];

        // Scale from degrees (radius 90) to pixel units (radius d)
        const xScaled = x.map(v => v * d / 90);
        const yScaled = y.map(v => v * d / 90);

        // Pixel resolution relative to canvas size
        const res = width.map(v => v / size);

        // Convert (x,y) arrays to a drawable patch and store
        collect.push(xy2patch(yScaled, xScaled, res, xy));
    }

    // Helper: numpy-like arange
    function range(start, stop, step = 1) {
        const out = [];
        if (step > 0) {
            for (let v = start; v < stop; v += step) out.push(v);
        } else {
            for (let v = start; v > stop; v += step) out.push(v);
        }
        return out;
    }

    return [['b', 'w'], collect];   // colour order + collection of polygons
}

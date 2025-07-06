// plotMT.js:

import { xy2patch } from "../math/geometry";
import { D2R, EPSILON } from "../math/trigHelpers";
 
/**
 * Generate patch objects and corresponding fill colors for a moment-tensor beachball.
 *
 * @param {PrincipalAxis} T             Tension principal axis.
 * @param {PrincipalAxis} N             Neutral/principal axis.
 * @param {PrincipalAxis} P             Pressure principal axis.
 * @param {number}        size          Diameter of the beachball in pixels.
 * @param {boolean}       plotZeroTrace If true, force zero isotropic component.
 * @param {number}        x0            X-coordinate offset of the center.
 * @param {number}        y0            Y-coordinate offset of the center.
 * @param {[number,number]} xy          Center position as [x, y].
 * @param {[number,number]|number} width Ellipse radii: single value or [width, height].
 * @returns {[string[], object[]]}       A tuple of:
 *   - colors: array of fill-color strings (e.g. `"b"`, `"w"`),  
 *   - patches: array of patch-definition objects for drawing.
 */
export function plotMT(T, N, P, size = 200, plot_zerotrace = true, x0 = 0, y0 = 0, xy = [0, 0], width = 200) {

    /* ---- check if one or two widths are specified (Circle / Ellipse) */
    /* from matplotlib import patches  ← irrelevant in JS               */
    try {
        if (!Array.isArray(width) || width.length !== 2) throw new TypeError();
    } catch (e) {
        width = [width, width];
    }

    /* ---- Python lists replaced with JS arrays of same length -------- */
    const collect = [];
    const colors = [];
    const res = width.map(value => value / size);
    let b = 1;
    let big_iso = 0;
    let j = 1;
    let j2 = 0;
    let j3 = 0;
    let n = 0;

    const azi = Array.from({ length: 3 }, () => new Array(2).fill(0));
    const x = new Array(400).fill(0);
    const y = new Array(400).fill(0);
    const x2 = new Array(400).fill(0);
    const y2 = new Array(400).fill(0);
    const x3 = new Array(400).fill(0);
    const y3 = new Array(400).fill(0);
    const xp1 = new Array(800).fill(0);
    const yp1 = new Array(800).fill(0);
    const xp2 = new Array(400).fill(0);
    const yp2 = new Array(400).fill(0);

    /* ---- a, p, v arrays -------------------------------------------- */
    const a = [T.strike, N.strike, P.strike];
    const p = [T.dip, N.dip, P.dip];
    const v = [T.val, N.val, P.val];

    let vi = (v[0] + v[1] + v[2]) / 3.0;
    for (let i = 0; i < 3; i++) {
        v[i] = v[i] - vi;
    }

    const radius_size = size * 0.5;

    /* ---- pure implosion-explosion check ----------------------------- */
    if (Math.abs(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) < EPSILON) {
        if (vi > 0.0) {
            const cir = { type: "ellipse", xy, width: width[0], height: width[1] };
            collect.push(cir);
            colors.push("b");
        }
        if (vi < 0.0) {
            const cir = { type: "ellipse", xy, width: width[0], height: width[1] };
            collect.push(cir);
            colors.push("w");
        }
        return [colors, collect];
    }

    /* ---- choose d / m ---------------------------------------------- */
    let d, m;
    if (Math.abs(v[0]) >= Math.abs(v[2])) { d = 0; m = 2; }
    else { d = 2; m = 0; }

    if (plot_zerotrace) vi = 0.0;

    const f = -v[1] / v[d];
    const iso = vi / v[d];

    // # Cliff Frohlich, Seismological Research letters,
    // # Vol 7, Number 1, January-February, 1996
    // # Unless the isotropic parameter lies in the range
    // # between -1 and 1 - f there will be no nodes whatsoever
    /* ---- Frohlich node existence test ------------------------------ */
    if (iso < -1) {
        collect.push({ type: "ellipse", xy, width: width[0], height: width[1] });
        colors.push("w");
        return [colors, collect];
    } else if (iso > 1 - f) {
        collect.push({ type: "ellipse", xy, width: width[0], height: width[1] });
        colors.push("b");
        return [colors, collect];
    }

    /* ---- lots of trig shorthands ----------------------------------- */
    const spd = Math.sin(p[d] * D2R);
    const cpd = Math.cos(p[d] * D2R);
    const spb = Math.sin(p[b] * D2R);
    const cpb = Math.cos(p[b] * D2R);
    const spm = Math.sin(p[m] * D2R);
    const cpm = Math.cos(p[m] * D2R);
    const sad = Math.sin(a[d] * D2R);
    const cad = Math.cos(a[d] * D2R);
    const sab = Math.sin(a[b] * D2R);
    const cab = Math.cos(a[b] * D2R);
    const sam = Math.sin(a[m] * D2R);
    const cam = Math.cos(a[m] * D2R);

    /* ---- main 0-to-359 loop ---------------------------------------- */
    let azp = 0.0;
    for (let i = 0; i < 360; i++) {
        const fir = i * D2R;
        const s2alphan = (2.0 + 2.0 * iso) /
            (3.0 + (1.0 - 2.0 * f) * Math.cos(2.0 * fir));

        if (s2alphan > 1.0) {
            big_iso += 1;
        } else {
            const alphan = Math.asin(Math.sqrt(s2alphan));
            const sfi = Math.sin(fir);
            const cfi = Math.cos(fir);
            const san = Math.sin(alphan);
            const can = Math.cos(alphan);

            const xz = can * spd + san * sfi * spb + san * cfi * spm;
            const xn = can * cpd * cad + san * sfi * cpb * cab +
                san * cfi * cpm * cam;
            const xe = can * cpd * sad + san * sfi * cpb * sab +
                san * cfi * cpm * sam;

            let takeoff, az;
            if (Math.abs(xn) < EPSILON && Math.abs(xe) < EPSILON) {
                takeoff = 0.0;
                az = 0.0;
            } else {
                az = Math.atan2(xe, xn);
                if (az < 0.0) az += Math.PI * 2.0;
                takeoff = Math.acos(
                    xz / Math.sqrt(xz * xz + xn * xn + xe * xe)
                );
            }

            if (takeoff > Math.PI / 2.0) {
                takeoff = Math.PI - takeoff;
                az += Math.PI;
                if (az > Math.PI * 2.0) az -= Math.PI * 2.0;
            }

            const r = Math.SQRT2 * Math.sin(takeoff / 2.0);
            const si = Math.sin(az);
            const co = Math.cos(az);

            if (i === 0) {
                azi[i][0] = az;
                x[i] = x0 + radius_size * r * si;
                y[i] = y0 + radius_size * r * co;
                azp = az;
            } else {
                if (Math.abs(Math.abs(az - azp) - Math.PI) < D2R * 10.0) {
                    azi[n][1] = azp;
                    n += 1;
                    azi[n][0] = az;
                }
                if (Math.abs(Math.abs(az - azp) - Math.PI * 2.0) < D2R * 2.0) {
                    if (azp < az) azi[n][0] += Math.PI * 2.0;
                    else azi[n][0] -= Math.PI * 2.0;
                }
                if (n === 0) {
                    x[j] = x0 + radius_size * r * si;
                    y[j] = y0 + radius_size * r * co;
                    j += 1;
                } else if (n === 1) {
                    x2[j2] = x0 + radius_size * r * si;
                    y2[j2] = y0 + radius_size * r * co;
                    j2 += 1;
                } else if (n === 2) {
                    x3[j3] = x0 + radius_size * r * si;
                    y3[j3] = y0 + radius_size * r * co;
                    j3 += 1;
                }
                azp = az;
            }
        }
    }
    azi[n][1] = azp;

    /* ---- color selections based on v[1] ----------------------------- */
    let rgb1, rgb2;
    if (v[1] < 0.0) {
        rgb1 = "b";
        rgb2 = "w";
    } else {
        rgb1 = "w";
        rgb2 = "b";
    }

    /* ---- background ellipse ---------------------------------------- */
    const cir = { type: "ellipse", xy, width: width[0], height: width[1] };
    collect.push(cir);
    colors.push(rgb2);

    /* ---- n == 0 case ------------------------------------------------ */
    if (n === 0) {
        collect.push(xy2patch(x.slice(0, 360), y.slice(0, 360), res, xy)); 
        colors.push(rgb1);
        return [colors, collect];
    }

    /* ---- n == 1 case ------------------------------------------------ */
    if (n === 1) {
        /* copy x,y into xp1,yp1 */
        let i;
        for (i = 0; i < j; i++) {
            xp1[i] = x[i];
            yp1[i] = y[i];
        }
        i -= 1;

        /* unwrap azi[0][0] w.r.t. azi[0][1] */
        if (azi[0][0] - azi[0][1] > Math.PI) {
            azi[0][0] -= Math.PI * 2.0;
        } else if (azi[0][1] - azi[0][0] > Math.PI) {
            azi[0][0] += Math.PI * 2.0;
        }

        /* add arc points between azi[0][1] and azi[0][0] */
        if (azi[0][0] < azi[0][1]) {
            let az = azi[0][1] - D2R;
            while (az > azi[0][0]) {
                const si = Math.sin(az);
                const co = Math.cos(az);
                xp1[i] = x0 + radius_size * si;
                yp1[i] = y0 + radius_size * co;
                i += 1;
                az -= D2R;
            }
        } else {
            let az = azi[0][1] + D2R;
            while (az < azi[0][0]) {
                const si = Math.sin(az);
                const co = Math.cos(az);
                xp1[i] = x0 + radius_size * si;
                yp1[i] = y0 + radius_size * co;
                i += 1;
                az += D2R;
            }
        }

        collect.push(xy2patch(xp1.slice(0, i), yp1.slice(0, i), res, xy));
        colors.push(rgb1);

        /* copy x2,y2 into xp2,yp2 */
        for (i = 0; i < j2; i++) {
            xp2[i] = x2[i];
            yp2[i] = y2[i];
        }
        i -= 1;

        /* unwrap azi[1][0] w.r.t azi[1][1] */
        if (azi[1][0] - azi[1][1] > Math.PI) {
            azi[1][0] -= Math.PI * 2.0;
        } else if (azi[1][1] - azi[1][0] > Math.PI) {
            azi[1][0] += Math.PI * 2.0;
        }

        /* arc for second node set */
        if (azi[1][0] < azi[1][1]) {
            let az = azi[1][1] - D2R;
            while (az > azi[1][0]) {
                const si = Math.sin(az);
                const co = Math.cos(az);
                xp2[i] = x0 + radius_size * si;
                i += 1;
                yp2[i] = y0 + radius_size * co;
                az -= D2R;
            }
        } else {
            let az = azi[1][1] + D2R;
            while (az < azi[1][0]) {
                const si = Math.sin(az);
                const co = Math.cos(az);
                xp2[i] = x0 + radius_size * si;
                i += 1;
                yp2[i] = y0 + radius_size * co;
                az += D2R;
            }
        }

        collect.push(xy2patch(xp2.slice(0, i), yp2.slice(0, i), res, xy));
        colors.push(rgb1);
        return [colors, collect];
    }

    /* ---- n == 2 case ------------------------------------------------ */
    if (n === 2) {
        /* first concatenate x3 → xp1 */
        let i = 0;
        for (i = 0; i < j3; i++) {
            xp1[i] = x3[i];
            yp1[i] = y3[i];
        }
        i -= 1;

        /* then x → xp1 */
        for (let ii = 0; ii < j; ii++) {
            xp1[i] = x[ii];
            i++;
            yp1[i] = y[ii];
        }

        /* big_iso shortcut --------------------------------------------- */
        if (big_iso) {
            let ii = j2 - 1;
            while (ii >= 0) {
                xp1[i] = x2[ii];
                i += 1;
                yp1[i] = y2[ii];
                ii -= 1;
            }
            collect.push(xy2patch(xp1.slice(0, i), yp1.slice(0, i), res, xy));
            colors.push(rgb1);
            return [colors, collect];
        }

        /* unwrap azi[2][0] wrt azi[0][1] */
        if (azi[2][0] - azi[0][1] > Math.PI) {
            azi[2][0] -= Math.PI * 2.0;
        } else if (azi[0][1] - azi[2][0] > Math.PI) {
            azi[2][0] += Math.PI * 2.0;
        }

        /* arc between azi[0][1] and azi[2][0] */
        if (azi[2][0] < azi[0][1]) {
            let az = azi[0][1] - D2R;
            while (az > azi[2][0]) {
                const si = Math.sin(az);
                const co = Math.cos(az);
                xp1[i] = x0 + radius_size * si;
                i += 1;
                yp1[i] = y0 + radius_size * co;
                az -= D2R;
            }
        } else {
            let az = azi[0][1] + D2R;
            while (az < azi[2][0]) {
                const si = Math.sin(az);
                const co = Math.cos(az);
                xp1[i] = x0 + radius_size * si;
                i += 1;
                yp1[i] = y0 + radius_size * co;
                az += D2R;
            }
        }

        collect.push(xy2patch(xp1.slice(0, i), yp1.slice(0, i), res, xy));
        colors.push(rgb1);

        /* copy x2 → xp2 */
        for (i = 0; i < j2; i++) {
            xp2[i] = x2[i];
            yp2[i] = y2[i];
        }
        i -= 1;

        /* unwrap azi[1][0] wrt azi[1][1] */
        if (azi[1][0] - azi[1][1] > Math.PI) {
            azi[1][0] -= Math.PI * 2.0;
        } else if (azi[1][1] - azi[1][0] > Math.PI) {
            azi[1][0] += Math.PI * 2.0;
        }

        /* arc for second patch */
        if (azi[1][0] < azi[1][1]) {
            let az = azi[1][1] - D2R;
            while (az > azi[1][0]) {
                const si = Math.sin(az);
                const co = Math.cos(az);
                xp2[i] = x0 + radius_size * si;
                i += 1;
                yp2[i] = y0 + radius_size * co;
                az -= D2R;
            }

        } else {
            let az = azi[1][1] + D2R;
            while (az < azi[1][0]) {
                const si = Math.sin(az);
                const co = Math.cos(az);
                xp2[i] = x0 + radius_size * si;
                i += 1;
                yp2[i] = y0 + radius_size * co;
                az += D2R;
            }
        }

        collect.push(xy2patch(xp2.slice(0, i), yp2.slice(0, i), res, xy));
        colors.push(rgb1);
        return [colors, collect];
    }

    /* ---- fallback (should never happen) ----------------------------- */
    return [colors, collect];

}

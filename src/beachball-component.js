// import necessary classes and helper
import { MomentTensor, NodalPlane } from './classes.js';
import { mt2plane, mt2axes, plotDC, plotMT, sdr2mt } from './beachball.js';
import { EPSILON } from './constants.js';

class BeachballComponent extends HTMLElement {
    static get observedAttributes() {
        return [
            'focal-mechanism',
            'width',
            'height',
            'tension-color',
            'bg-color',
            'edge-color',
            'line-width',
            'size'
        ];
    }

    constructor() {
        super();
        const shadow = this.attachShadow({ mode: 'open' });

        // styles for host & canvas
        const style = document.createElement('style');
        style.textContent = `
      :host {
        display: inline-block;
        position: relative;
        width: var(--bb-width, 200px);
        height: var(--bb-height, 200px);
      }
      canvas {
        width: 100%;
        height: 100%;
        display: block;
      }
    `;
        shadow.appendChild(style);

        // only the canvas
        this.canvas = document.createElement('canvas');
        shadow.appendChild(this.canvas);
    }

    attributeChangedCallback(name, oldVal, newVal) {
        switch (name) {
            case 'focal-mechanism':
                try {
                    this.focalMechanism = JSON.parse(newVal);
                } catch {
                    console.warn('Invalid JSON for focal-mechanism');
                }
                break;
            case 'width': {
                const w = parseInt(newVal) || 200;
                this.canvas.width = w;
                if (!this.hasAttribute('height')) {
                    this.canvas.height = w;
                    this.style.setProperty('--bb-height', `${w}px`);
                }
                this.style.setProperty('--bb-width', `${w}px`);
                break;
            }
            case 'height': {
                const h = parseInt(newVal);
                this.canvas.height = h;
                this.style.setProperty('--bb-height', `${h}px`);
                break;
            }
            case 'tension-color':
                this.tensionColor = newVal;
                break;
            case 'bg-color':
                this.backgroundColor = newVal;
                break;
            case 'edge-color':
                this.edgeColor = newVal;
                break;
            case 'line-width':
                this.lineWidth = parseFloat(newVal);
                break;
            case 'size':
                this.size = Math.max(100, parseInt(newVal));
                break;
        }
        this._draw();
    }

    set focalMechanism(fm) {
        let mt = null, np1 = null;
        if (Array.isArray(fm) && fm.length === 6) {
            mt = new MomentTensor(...fm, 0);
            np1 = mt2plane(mt);
        } else if (Array.isArray(fm) && fm.length === 3) {
            np1 = new NodalPlane(...fm);
        } else {
            throw new TypeError("Wrong input value for 'fm'.");
        }
        this._mt = mt;
        this._np1 = np1;
        this._draw();
    }
    get focalMechanism() {
        return this._mt ? this._mt : this._np1;
    }

    connectedCallback() {
        // defaults
        const w = this.hasAttribute('width') ? +this.getAttribute('width') : 200;
        this.canvas.width = w;
        this.style.setProperty('--bb-width', `${w}px`);

        const h = this.hasAttribute('height') ? +this.getAttribute('height') : w;
        this.canvas.height = h;
        this.style.setProperty('--bb-height', `${h}px`);

        this.tensionColor = this.getAttribute('tension-color') || '#87CEEB';
        this.backgroundColor = this.getAttribute('bg-color') || '#FFFFFF';
        this.edgeColor = this.getAttribute('edge-color') || '#000000';
        this.lineWidth = +this.getAttribute('line-width') || 2;
        this.size = Math.max(100, +this.getAttribute('size') || 100);

        this._draw();
    }

    _drawPatches(ctx, colors, patches, fill = true, drawLine = false) {
        ctx.save();

        // --- flip y-axis once for ALL patches ---
        // Adjust coordinate system: Matplotlib uses a bottom-left origin with y-axis up,
        // whereas Canvas uses a top-left origin with y-axis down. This translates the
        // origin to the bottom and flips the y-axis so patches render with correct symmetry.
        ctx.translate(0, ctx.canvas.height);
        ctx.scale(1, -1);

        ctx.lineWidth = this.lineWidth;
        ctx.strokeStyle = this.edgeColor;

        patches.forEach((patch, i) => {
            ctx.fillStyle = colors[i];
            ctx.beginPath();
            if (patch.type === 'ellipse') {
                ctx.ellipse(
                    patch.xy[0], patch.xy[1],
                    patch.width / 2, patch.height / 2,
                    0, 0, 2 * Math.PI
                );

            } else {
                patch.vertices.forEach(([x, y], idx) =>
                    idx ? ctx.lineTo(x, y) : ctx.moveTo(x, y)
                );
            }
            ctx.closePath();
            if (fill)
                ctx.fill();
            if (drawLine)
                ctx.stroke();
        });
        ctx.restore();
    }


    /**
     * Project an axis onto a beachball diagram.
     *
     * @param {number} strike – strike angle in degrees
     * @param {number} dip    – dip angle in degrees
     * @param {number} [R=90] – reference radius (pixels or arbitrary units)
     * @param {number} [margin=20] – extra margin added to the projected radius
     * @returns {[number, number]} [x, y] coordinates in the plotting plane
     */
    _axisOnBeachball(strike, dip, R = 90, margin = 20) {
        const az = strike * Math.PI / 180;              // radians
        const r = (90 - dip) * R / 90 + margin;        // projected radius
        const y = Math.cos(az) * r;
        const x = Math.sin(az) * r;
        return [x, y];                                  // same order as Python
    }


    /**
    * Return an array of fill colors based on the raw `colors` array:
    * if colors[i] === 'b' use tensionColor, else backgroundColor.
    */
    _getFillColors(colors, patches) {
        return patches.map((_, i) =>
            colors[i] === 'b' ? this.tensionColor : this.backgroundColor
        );
    }

    _draw() {
        if (!this._np1) return;
        const ctx = this.canvas.getContext('2d');

        // fill background
        ctx.fillStyle = this.backgroundColor;
        ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);

        // prepare border style
        ctx.lineWidth = this.lineWidth;   // e.g. 2
        ctx.strokeStyle = this.edgeColor;    // e.g. 'black'

        // compute center and radius (or radii)
        const cx = this.canvas.width / 2;
        const cy = this.canvas.height / 2;

        // // draw the main border (circle or ellipse)
        // ctx.beginPath();
        // ctx.ellipse(cx, cy, cx - ctx.lineWidth, cy - ctx.lineWidth, 0, 0, 2 * Math.PI);
        // ctx.stroke();

        // now draw internal patches
        const xy = [cx, cy];
        const diameter = [this.canvas.width - 2 * ctx.lineWidth, this.canvas.height - 2 * ctx.lineWidth]; //r * 2;
        let colors, patches;

        if (this._mt) {
            const [t, n, p] = mt2axes(this._mt.normalized);

            // radius of the plotting circle in canvas units:
            const R = Math.min(this.canvas.width, this.canvas.height) / 2 - this.lineWidth;
            // same margin you used in Python:
            const MARGIN = 0; //18

            // 1) project to 2-D coordinates (matplotlib frame)
            const [xT, yT] = this._axisOnBeachball(t.strike, t.dip, R, MARGIN);
            const [xP, yP] = this._axisOnBeachball(p.strike, p.dip, R, MARGIN);

            if (Math.abs(n.val) < EPSILON && Math.abs(t.val + p.val) < EPSILON) {
                // pure DC
                [colors, patches] = plotDC(this._np1, this.size, xy, diameter);
                const fillColors = this._getFillColors(colors, patches);
                this._drawPatches(ctx, fillColors, patches, true, true);

            } else {
                // MT + DC overlay
                [colors, patches] = plotMT(t, n, p, this.size, true, 0, 0, xy, diameter);
                const fillMT = this._getFillColors(colors, patches);
                this._drawPatches(ctx, fillMT, patches, true, false);

                [colors, patches] = plotDC(this._np1, this.size, xy, diameter);
                const fillDC = this._getFillColors(colors, patches);
                this._drawPatches(ctx, fillDC, patches, false, true);
            }

            // 2) flip to canvas coordinates (because we use y-down on canvas)
            const canvasXT = cx + xT;
            const canvasYT = cy - yT;   // minus sign flips y
            const canvasXP = cx + xP;
            const canvasYP = cy - yP;

            // 3) draw bold centred text
            ctx.save();
            ctx.font = 'bold 18px sans-serif';
            ctx.fillStyle = 'black';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            // “z-order” on canvas == call order, so draw last
            ctx.fillText('T', canvasXT, canvasYT);
            ctx.fillText('P', canvasXP, canvasYP);

            ctx.restore();


            console.log(">>>>>>>>>>>>>>>>>>>>>>>>>>>>");
            let mt2 = sdr2mt(240, 46, -119);
            console.log(mt2);




        } else {
            // fallback to DC only
            [colors, patches] = plotDC(this._np1, this.size, xy, diameter);
            const fillColors = this._getFillColors(colors, patches);
            this._drawPatches(ctx, fillColors, patches, true, true);
        }
    }


}

customElements.define('beachball-component', BeachballComponent);

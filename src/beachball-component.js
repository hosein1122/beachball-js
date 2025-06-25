// beachball-component.js
// Web Component for rendering beachball focal mechanisms on a HTML canvas.

import { beach } from "./beachball.js";
import { patchToPath2D } from "./canvasUtils.js";

class BeachballViewer extends HTMLElement {
    static get observedAttributes() {
        return [
            'fm',           // JSON array: [strike,dip,rake] or [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp]
            'facecolor',    // color for tension quadrants
            'bgcolor',      // color for pressure quadrants
            'edgecolor',    // outline color
            'alpha',        // fill opacity
            'radius',       // radius in px
            'radius-y'      // optional y-axis radius
        ];
    }

    constructor() {
        super();
        this.attachShadow({ mode: 'open' });
        this.canvas = document.createElement('canvas');
        this.shadowRoot.appendChild(this.canvas);
    }

    connectedCallback() {
        this.updateSize();
        window.addEventListener('resize', () => { this.updateSize(); this.render(); });
        this.render();
    }

    attributeChangedCallback() {
        this.render();
    }

    updateSize() {
        const rect = this.getBoundingClientRect();
        this.canvas.width = rect.width;
        this.canvas.height = rect.height;
    }

    getProps() {
        let fmAttr = this.getAttribute('fm');
        let fm = [];
        try { fm = JSON.parse(fmAttr); } catch { };

        const facecolor = this.getAttribute('facecolor') || 'b';
        const bgcolor = this.getAttribute('bgcolor') || 'w';
        const edgecolor = this.getAttribute('edgecolor') || 'k';
        const alpha = parseFloat(this.getAttribute('alpha')) || 1.0;
        const radius = parseFloat(this.getAttribute('radius')) || Math.min(this.canvas.width, this.canvas.height) / 2;
        const radiusY = parseFloat(this.getAttribute('radius-y')) || radius;

        return { fm, facecolor, bgcolor, edgecolor, alpha, radius, radiusY };
    }

    render() {
        const ctx = this.canvas.getContext('2d');
        ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        ctx.save();
        // center origin
        ctx.translate(this.canvas.width / 2, this.canvas.height / 2);

        const { fm, facecolor, bgcolor, edgecolor, alpha, radius, radiusY } = this.getProps();

        // Determine whether 6-component moment-tensor or 3-component nodal plane
        const isMT = Array.isArray(fm) && fm.length === 6;
        const opts = { facecolor, bgcolor, edgecolor, alpha, radius, radiusY };

        let patches;
        if (isMT) {
            patches = beach(fm, opts);
        } else {
            // 3-component: only double-couple
            opts.nofill = false;
            patches = beach(fm, opts);
        }

        // Helper to draw a patch
        const draw = patch => {
            const path = patchToPath2D(patch);
            ctx.globalAlpha = patch.alpha;
            ctx.fillStyle = patch.facecolor;
            ctx.strokeStyle = patch.edgecolor;
            if (patch.vertices.length) ctx.fill(path);
            ctx.stroke(path);
        };

        draw(patches.pressure);
        draw(patches.tension);

        ctx.restore();
    }
}

customElements.define('beachball-viewer', BeachballViewer);

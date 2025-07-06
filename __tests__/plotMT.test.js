// __tests__/plotMT.test.js
import { PrincipalAxis } from "../src/classes.js";
import { plotMT } from "../src/beachball/plotMT";

describe("plotMT (pure‐isotropic cases)", () => {
    test("pure explosion (positive isotropic) → one blue ellipse", () => {
        // T=N=P all equal ⇒ pure isotropic positive
        const T = new PrincipalAxis(1.0, 0, 0);
        const N = new PrincipalAxis(1.0, 0, 0);
        const P = new PrincipalAxis(1.0, 0, 0);

        const [colors, patches] = plotMT(T, N, P, /*size*/100);
        expect(colors).toEqual(["b"]);
        expect(patches).toHaveLength(1);

        const patch = patches[0];
        expect(patch.type).toBe("ellipse");
        expect(patch.xy).toEqual([0, 0]);
        expect(patch.width).toBe(200);
        expect(patch.height).toBe(200);
    });

    test("pure implosion (negative isotropic) → one white ellipse", () => {
        // T=N=P all equal negative ⇒ pure isotropic negative
        const T = new PrincipalAxis(-2.0, 0, 0);
        const N = new PrincipalAxis(-2.0, 0, 0);
        const P = new PrincipalAxis(-2.0, 0, 0);

        const [colors, patches] = plotMT(T, N, P, /*size*/50);
        expect(colors).toEqual(["w"]);
        expect(patches).toHaveLength(1);

        const patch = patches[0];
        expect(patch.type).toBe("ellipse");
        expect(patch.xy).toEqual([0, 0]);
        expect(patch.width).toBe(200);
        expect(patch.height).toBe(200);
    });
});

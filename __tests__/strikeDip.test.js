import { strikeDip } from "../src/math/geometry.js";

describe("strikeDip", () => {
    test("horizontal plane (normal up) → strike 270°, dip 0°", () => {
        const [strike, dip] = strikeDip(0, 0, 1);
        expect(strike).toBeCloseTo(270, 6);
        expect(dip).toBeCloseTo(0, 6);
    });

    test("vertical plane facing north (n=1,e=0,u=0) → strike 270°, dip 90°", () => {
        const [strike, dip] = strikeDip(1, 0, 0);
        expect(strike).toBeCloseTo(270, 6);
        expect(dip).toBeCloseTo(90, 6);
    });

    test("vertical plane facing east (n=0,e=1,u=0) → strike 0°, dip 90°", () => {
        const [strike, dip] = strikeDip(0, 1, 0);
        expect(strike).toBeCloseTo(0, 6);
        expect(dip).toBeCloseTo(90, 6);
    });

    test("downward normal flipped to upward (0,0,-1) → same as (0,0,1)", () => {
        const up = strikeDip(0, 0, 1);
        const down = strikeDip(0, 0, -1);
        expect(down[0]).toBeCloseTo(up[0], 6);
        expect(down[1]).toBeCloseTo(up[1], 6);
    });

    test("normal (1,1,1) → strike ≈315°, dip ≈54.7356°", () => {
        const norm = Math.sqrt(3);
        const [strike, dip] = strikeDip(1 / norm, 1 / norm, 1 / norm);
        expect(strike).toBeCloseTo(315, 3);
        expect(dip).toBeCloseTo(54.7356, 4);
    });
});

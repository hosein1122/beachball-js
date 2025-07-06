// __tests__/auxPlane.test.js
import { auxPlane } from "../src/math/geometry.js";

describe("auxPlane", () => {
    test("horizontal plane (s1=0, d1=0, r1=0) → [270°, 90°, –90°]", () => {
        const [strike2, dip2, rake2] = auxPlane(0, 0, 0);
        expect(strike2).toBeCloseTo(270, 6);
        expect(dip2).toBeCloseTo(90, 6);
        expect(rake2).toBeCloseTo(-90, 6);
    });

    test("vertical N–S plane (s1=0, d1=90, r1=0) → [270°, 90°, –180°]", () => {
        const [strike2, dip2, rake2] = auxPlane(0, 90, 0);
        expect(strike2).toBeCloseTo(270, 6);
        expect(dip2).toBeCloseTo(90, 6);
        expect(rake2).toBeCloseTo(-180, 6);
    });

    test(" (s1=45, d1=30, r1=45) → [274.106605°, 69.295189°, 112.207654°]", () => {
        const [strike2, dip2, rake2] = auxPlane(45, 30, 45);
        expect(strike2).toBeCloseTo(274.106605, 6);
        expect(dip2).toBeCloseTo(69.295189, 6);
        expect(rake2).toBeCloseTo(112.207654, 6);
    });


    // you can add more cases here, e.g. based on known beachball examples
});

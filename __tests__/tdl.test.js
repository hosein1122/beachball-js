// __tests__/tdl.test.js
import { tdl } from "../src/math/geometry.js";

describe("tdl", () => {
    test("flat‐horizontal branch: an=[1,0,0], bn=[0,0,1] → [ft,fd,fl] = [270,90,−90]", () => {
        const [ft, fd, fl] = tdl([1, 0, 0], [0, 0, 1]);
        expect(ft).toBeCloseTo(270, 6);
        expect(fd).toBeCloseTo(90, 6);
        expect(fl).toBeCloseTo(-90, 6);
    });

    test("an=[0, 0, 1], [1, 0, 0] → [ft,fd,fl] = [0,180,0]", () => {
        const [ft, fd, fl] = tdl([0, 0, 1], [1, 0, 0]);
        expect(ft).toBeCloseTo(0, 6);
        expect(fd).toBeCloseTo(180, 6);
        expect(fl).toBeCloseTo(0, 6);
    });

});

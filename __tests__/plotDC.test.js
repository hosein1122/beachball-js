// __tests__/plotDC.test.js
import { plotDC } from "../src/beachball/plotDC.js";

describe("plotDC", () => {
    test("horizontal plane (strike=0,dip=0,rake=0) → basic structure", () => {
        const np1 = { strike: 0, dip: 0, rake: 0 };
        const size = 200;
        const centre = [0, 0];
        const width = 200; // circle

        const [colours, patches] = plotDC(np1, size, centre, width);

        // 1) رنگ‌ها و تعداد پچ‌ها
        expect(colours).toEqual(["b", "w"]);
        expect(patches).toHaveLength(2);

        // 2) هر پچ باید شامل vertices, res, center باشد
        patches.forEach(patch => {
            expect(Array.isArray(patch.vertices)).toBe(true);
            expect(patch.vertices.length).toBeGreaterThan(0);
            patch.vertices.forEach(pt => {
                expect(Array.isArray(pt)).toBe(true);
                expect(pt).toHaveLength(2);
                expect(typeof pt[0]).toBe("number");
                expect(typeof pt[1]).toBe("number");
            });

            expect(Array.isArray(patch.res)).toBe(true);
            expect(patch.res).toHaveLength(2);
            expect(patch.res).toEqual([width / size, width / size]);

            expect(Array.isArray(patch.center)).toBe(true);
            expect(patch.center).toEqual(centre);
        });
    });

    test("ellipse width array scales axes differently", () => {
        const np1 = { strike: 45, dip: 30, rake: 60 };
        const size = 100;
        const centre = [10, -20];
        const width = [50, 80]; // ellipse: rx=50, ry=80

        const [_, patches] = plotDC(np1, size, centre, width);
        patches.forEach(patch => {
            // res should be [50/100, 80/100]
            expect(patch.res).toEqual([0.5, 0.8]);
            // centre unchanged
            expect(patch.center).toEqual(centre);
        });
    });

    test("invalid width length throws", () => {
        const np1 = { strike: 0, dip: 0, rake: 0 };
        expect(() => {
            plotDC(np1, 100, [0, 0], [10]);
        }).toThrow("width must contain exactly two elements");
    });
});

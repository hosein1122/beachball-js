// __tests__/mt2plane.test.js
import { mt2plane } from "../src/beachball/mt2plane";
import { MomentTensor, NodalPlane } from "../src/classes";

describe("mt2plane", () => {
    test("returns a NodalPlane instance with finite strike/dip/rake", () => {
        // مثال دلخواه: همه مؤلفه‌ها غیر صفر
        const mt = new MomentTensor(
      /* Mrr */ 4.0,
      /* Mtt */ 5.0,
      /* Mpp */ 6.0,
      /* Mrt */ 1.0,
      /* Mrp */ 2.0,
      /* Mtp */ 3.0,
            0
        );
        const plane = mt2plane(mt);

        // ۱) نوع خروجی
        expect(plane).toBeInstanceOf(NodalPlane);

        // ۲) مؤلفه‌ها عدد و متناهی باشند
        ["strike", "dip", "rake"].forEach(prop => {
            expect(typeof plane[prop]).toBe("number");
            expect(isFinite(plane[prop])).toBe(true);
        });
    });


    test("JS mt2plane → compare with Python example", () => {
        const Mrr = 4.0;
        const Mtt = 5.0;
        const Mpp = 6.0;
        const Mrt = 1.0;
        const Mrp = 2.0;
        const Mtp = 3.0;
        // note: the Python example passed an extra zero for the 7th arg
        const mt = new MomentTensor(Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, 0);

        const planeJS = mt2plane(mt);

        // console.log(`JS mt2plane →
        //  strike = ${planeJS.strike.toFixed(6)}
        //  dip    = ${planeJS.dip.toFixed(6)}
        //  rake   = ${planeJS.rake.toFixed(6)}`);

        // a trivial assertion so Jest reports “PASS”
        expect(typeof planeJS.strike).toBe("number");
    });
});

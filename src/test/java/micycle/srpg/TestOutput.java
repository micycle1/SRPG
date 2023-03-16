package micycle.srpg;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.random.RandomGenerator;

import org.junit.jupiter.api.Test;

/**
 * Test output of Java port vs original C implementation.
 * 
 * @author Michael Carleton
 *
 */
class TestOutput {

	@Test
	void testSimpleAligned() {
		int nX = 17;
		int nY = 15;
		double markPercent = 0.15;
		boolean holes = false;
		boolean aligned = false;
		boolean perturb = false;
		int smooth = 0;
		int hierarchy = 3;
		boolean diagonal = false;
		runCase(nX, nY, markPercent, holes, aligned, perturb, smooth, hierarchy, diagonal);
	}

	@Test
	void testPerturbedOrthogonal() {
		int nX = 20;
		int nY = 20;
		double markPercent = 0.25;
		boolean holes = false;
		boolean aligned = false;
		boolean perturb = true;
		int smooth = 0;
		int hierarchy = 2;
		boolean diagonal = false;
		runCase(nX, nY, markPercent, holes, aligned, perturb, smooth, hierarchy, diagonal);
	}
	
	@Test
	void testOrthogonal() {
		int nX = 12;
		int nY = 11;
		double markPercent = 0.33;
		boolean holes = false;
		boolean aligned = false;
		boolean perturb = false;
		int smooth = 0;
		int hierarchy = 4;
		boolean diagonal = false;
		runCase(nX, nY, markPercent, holes, aligned, perturb, smooth, hierarchy, diagonal);
	}
	
	@Test
	void testAlignedOrthogonal() {
		int nX = 10;
		int nY = 11;
		double markPercent = 0.45;
		boolean holes = false;
		boolean aligned = true;
		boolean perturb = false;
		int smooth = 0;
		int hierarchy = 3;
		boolean diagonal = false;
		runCase(nX, nY, markPercent, holes, aligned, perturb, smooth, hierarchy, diagonal);
	}
	
	@Test
	void testSmooth() {
		int nX = 10;
		int nY = 10;
		double markPercent = 0.05;
		boolean holes = false;
		boolean aligned = false;
		boolean perturb = true;
		int smooth = 3;
		int hierarchy = 3;
		boolean diagonal = false;
		runCase(nX, nY, markPercent, holes, aligned, perturb, smooth, hierarchy, diagonal);
	}
	
	@Test
	void testBigGrid() {
		int nX = 100;
		int nY = 100;
		double markPercent = 0.5;
		boolean holes = true; // NOTE
		boolean aligned = false;
		boolean perturb = true;
		int smooth = 2;
		int hierarchy = 0;
		boolean diagonal = true;
		runCase(nX, nY, markPercent, holes, aligned, perturb, smooth, hierarchy, diagonal);
	}

	void runCase(int Nx, int Ny, double mark_percent, boolean holes, boolean aligned, boolean perturb, int smoothRounds,
			int hierarchy, boolean diagonal) {
		String testCaseName = "%s,%s,%s,%s,%s,%s,%s,%s,%s".formatted(Nx, Ny, mark_percent, holes ? "t" : "f", aligned ? "t" : "f",
				perturb ? "t" : "f", smoothRounds, hierarchy, diagonal ? "t" : "f");

		var target = parseTestFile(testCaseName);
		SRPolygonGenerator srpg = new SRPolygonGenerator(Nx, Ny, mark_percent, holes, aligned, perturb, smoothRounds, hierarchy,
				diagonal, new XSRandom(1337));
		var actual = srpg.getPolygon();
		compare(target, actual);
	}

	void compare(List<List<double[]>> target, List<List<double[]>> actual) {
		assertEquals(target.size(), actual.size(), "Ring mismatch"); // compare number of rings
		for (int i = 0; i < target.size(); i++) {
			List<double[]> targetRing = target.get(i);
			List<double[]> actualRing = actual.get(i);
			assertEquals(targetRing.size(), actualRing.size(), "Ring #%s vertex number mismatch.".formatted(i)); // compare number of points in the rings
			for (int j = 0; j < targetRing.size(); j++) {
				double[] targetPoint = targetRing.get(j);
				double[] actualPoint = actualRing.get(j);
				assertArrayEquals(targetPoint, actualPoint, 1e-6, "Mismatch of coordinate #%s on ring #%s".formatted(j, i));
			}
		}
	}

	/**
	 * @return list polygon vertex rings provided by a given test file
	 */
	private static List<List<double[]>> parseTestFile(String name) {
		name+=".txt";
		Path filePath = Paths.get("src", "test", "resources", name);
		List<String> lines;
		try {
			lines = Files.readAllLines(filePath);
		} catch (IOException e) {
			assertTrue(false, "Could not find file: %s".formatted(name));
			return null;
		}
		List<List<double[]>> rings = new ArrayList<>();
		List<double[]> currentRing = new ArrayList<>();
		for (String line : lines) {
			line = line.trim();
			if (line.isEmpty()) {
				if (!currentRing.isEmpty()) {
					rings.add(currentRing);
				}
				currentRing = new ArrayList<>();
			} else {
				String[] parts = line.split(" ");
				double[] vertex = new double[2];
				vertex[0] = Double.parseDouble(parts[0]);
				vertex[1] = Double.parseDouble(parts[1]);
				currentRing.add(vertex);
				if (line == lines.get(lines.size() - 1)) {
					rings.add(currentRing);
				}
				continue;
			}
		}
		return rings;
	}

	/**
	 * Simple XorShift RandomGenerator. Used to ensures Java and C versions behave
	 * identically.
	 */
	private static class XSRandom implements RandomGenerator {
		private long seed;

		XSRandom(long seed) {
			this.seed = seed;
		}

		@Override
		public long nextLong() {
			long x = seed;
			x ^= (x << 21);
			x ^= (x >> 35);
			x ^= (x << 4);
			seed = x;
			return x;
		}

		@Override
		public int nextInt() {
			return (int) (Math.abs(nextLong()) % 2147483647);
		}
	}

}

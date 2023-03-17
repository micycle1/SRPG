package micycle.srpg;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.random.RandomGenerator;

/**
 * SRPG - Super Random Polygon Generator
 * <p>
 * SRPG generates simply-connected and multiply-connected polygons by means of a
 * regular grid that consists of square cells. Given two integer values, a and
 * b, SRPolygonGenerator generates a grid of size a times b.
 * <P>
 * By default SRPolygonGenerator then generates orthogonal polygons on this
 * grid. An additional parameter p, between zero and one, leads to a smaller or
 * larger number of vertices in the produced polygon. SRPolygonGenerator is able
 * to produce octagonal polygons by cutting off corners with ±45° diagonals
 * during the construction. Cutting corners repeatedly, without the diagonal
 * restriction, yields an approximation of a smooth free-form curve.
 * Additionally, SRPolygonGenerator can apply perturbations in order to generate
 * polygons with axes-parallel edges whose vertices do not lie on a grid, or to
 * generate polygons whose edges (in general) are not parallel to the coordinate
 * axes.
 * <p>
 * See
 * https://www1.pub.informatik.uni-wuerzburg.de/eurocg2020/data/uploads/papers/eurocg20_paper_75.pdf
 * for more info.
 * 
 * @author Michael Carleton (Java implementation)
 * @author Martin Held (algorithm)
 *
 */
public class SRPolygonGenerator {

	// https://github.com/cgalab/genpoly-micycle.srpg/blob/master/micycle.srpg.c
	// https://sbgdb.cs.sbg.ac.at/classes/polygons/
	// https://www1.pub.informatik.uni-wuerzburg.de/eurocg2020/data/uploads/papers/eurocg20_paper_75.pdf

	private static final int NIL = -1;

	private boolean[][] full;
	private boolean[][] keep;

	private boolean[][] top;
	private boolean[][] bot;
	private boolean[][] lft;
	private boolean[][] rgt;

	private int numCandidates = 0;
	private int numKeepCandidates = 0;

	private EdgeNode[][] edges;

	private int numVertices = 0;
	private int maxNumVertices = 0;

	private int nX;
	private int nY;
	private final boolean perturb;
	private boolean aligned;
	private final boolean diagonal;
	private int hierarchy = 0;
	private int smooth = 0;
	private double markPercent = 0.5;
	private boolean holes = false;

	private final RandomGenerator rand;

	private List<Coord> pnts = new ArrayList<>();
	private final List<Integer> candidates = new ArrayList<>();
	private final List<Integer> keepCandidates = new ArrayList<>();
	private VertexNode[] vertices = new VertexNode[0];

	private int maxKeep, maxCells;
	private int numKeep = 0;
	private int numCells = 0;

	/**
	 * Generates a random polygon based on a grid with Nx times Ny quadratic cells.
	 * The number of vertices of the polygon generated is random, but it does depend
	 * on Nx, Ny and the percentage markPercent: The larger Nx and Ny, the more
	 * vertices the polygon tends to have for a given markPercent.
	 * <p>
	 * Visually pleasing "random" polygons can be achieved by selecting fairly small
	 * values for markPercent, e.g., markPercent:=0.1 or even markPercent:=0.01.
	 * (However, a small value of markPercent will also reduce the number of
	 * vertices of the polygon.)
	 * 
	 * @param Nx           The number of cells in the X direction of the grid.
	 * @param Ny           The number of cells in the Y direction of the grid.
	 * @param markPercent  The percentage of vertices marked on the grid. The larger
	 *                     the percentage, the more vertices the polygon tends to
	 *                     have.
	 * @param holes        If true, generates a multiply-connected polygonal area.
	 * @param aligned      If true, all vertices will lie on grid points.
	 * @param perturb      If true, the vertices are moved away from the grid points
	 *                     and (most) polygon edges will not be parallel to the
	 *                     coordinate axes.
	 * @param smoothRounds The number of rounds of corner cutting (Chaikin) to apply
	 *                     to the polygon generated. A small positive integer value
	 *                     is recommended. A value of 3 is probably sufficient.
	 * @param hierarchy    The number of rounds of recursive refinement to apply to
	 *                     the polygon generated. A small positive integer value is
	 *                     recommended. This is akin to increasing the depth of
	 *                     fractal curve.
	 * @param diagonal     If true, cuts off some corners by line segments with
	 *                     inclination +/-1, thus generating an octagonal polygon.
	 */
	public SRPolygonGenerator(int Nx, int Ny, double mark_percent, final boolean holes, boolean aligned, boolean perturb, int smoothRounds,
			int hierarchy, final boolean diagonal, final RandomGenerator rand) {

		if (aligned) {
			perturb = false;
		} else if (perturb) {
			aligned = false;
		}
		if (hierarchy < 0) {
			hierarchy = 0;
		}
		if (smoothRounds < 0) {
			smoothRounds = 0;
		}
		if (Nx < 3) {
			Nx = 3;
		}
		if (Ny < 3) {
			Ny = 3;
		}
		if (mark_percent < 0.0001) {
			mark_percent = 0.0001;
		} else if (mark_percent > 0.5) {
			mark_percent = 0.5;
		}
		if (diagonal && !perturb) {
			this.aligned = true;
		}

		this.nX = Nx;
		this.nY = Ny;
		this.markPercent = mark_percent;
		this.holes = holes;
		this.aligned = aligned;
		this.perturb = perturb;
		this.smooth = smoothRounds;
		this.hierarchy = hierarchy;
		this.diagonal = diagonal;

		this.rand = rand;
	}

	/**
	 * 
	 * @return a random polygon generated based on the input parameters
	 */
	public List<List<double[]>> getPolygon() {
		return compute();
	}

	// generates the polygon
	private List<List<double[]>> compute() {
		// 1. allocate the grid
		int N = nX * nY;
		initializeGrid();

		// 2. select cells to be kept
		final double keepPercent = (1.0 - markPercent) * 0.95;
		maxKeep = (int) (N * keepPercent);
		maxCells = (int) (N * markPercent);

		int m = nY / 10;
		int k;
		for (int i = 0; i < nX; ++i) {
			k = uniformRandom(m);
			for (int j = 0; j <= k; ++j) {
				keep(i, j, false);
			}
			k = uniformRandom(m);
			for (int j = nY - 1; j >= nY - k - 1; --j) {
				keep(i, j, false);
			}
		}
		for (int j = 0; j < nY; ++j) {
			k = uniformRandom(m);
			for (int i = 0; i <= k; ++i) {
				keep(i, j, false);
			}
			k = uniformRandom(m);
			for (int i = nX - 1; i >= nX - k - 1; --i) {
				keep(i, j, false);
			}
		}

		int i, j, n;

		if (numKeep < maxKeep) {
			m = (maxKeep - numKeep) / 4;
			n = 0;
			for (n = 0; n < m; ++n) {
				k = uniformRandom(N);
				i = k / nY;
				j = k - i * nY;
				keep(i, j, true);
			}
		}

		i = uniformRandom(nX);
		j = uniformRandom(nY);

		if (numKeep < maxKeep) {
			for (int ii = i - nX / 30; ii < (i + nX / 30); ++ii) {
				for (int jj = j - nY / 30; jj < (j + nY / 30); ++jj) {
					keep(ii, jj, true);
					if (numKeep >= maxKeep) {
						break;
					}
				}
				if (numKeep >= maxKeep) {
					break;
				}
			}
		}

		// 3. select seed cell(s)
		m = 7 * nX / 9;
		n = 7 * nY / 9;
		final int mm = nX / 9;
		final int nn = nY / 9;
		do {
			i = uniformRandom(m);
			j = uniformRandom(n);
			i += mm;
			j += nn;
		} while (!isPossible(i, j));
		mark(i, j);

		// 4. randomly add cells to already selected cells
		selectRandomCells();

		// 5. use current polygon as "seed" for a refinement
		int hierarchyCntr = 1;
		boolean[][] old_full;
		while (hierarchyCntr <= hierarchy) {
			// 5.1 reset grid data
			old_full = full;
			final int N_x_old = nX;
			final int N_y_old = nY;
			nX *= 3;
			nY *= 3;
			N = nX * nY;

			initializeGrid();
			numCandidates = 0;
			numKeepCandidates = 0;
			maxKeep = (int) (N * keepPercent);
			maxCells = (int) (N * markPercent);
			numKeep = 0;
			numCells = 0; // =0

			// 5.2 copy data from coarse grid to fine grid
			for (i = 0; i < N_x_old; ++i) {
				for (j = 0; j < N_y_old; ++j) {
					final int ii = i * 3;
					final int jj = j * 3;
					if (isOldFull(i, j, N_x_old, N_y_old, old_full)) {
						mark(ii + 1, jj + 1);
						if (isOldFull(i - 1, j, N_x_old, N_y_old, old_full) && isOldFull(i, j - 1, N_x_old, N_y_old, old_full)) {
							mark(ii, jj);
						}
						if (isOldFull(i - 1, j, N_x_old, N_y_old, old_full) && isOldFull(i, j + 1, N_x_old, N_y_old, old_full)) {
							mark(ii, jj + 2);
						}
						if (isOldFull(i + 1, j, N_x_old, N_y_old, old_full) && isOldFull(i, j + 1, N_x_old, N_y_old, old_full)) {
							mark(ii + 2, jj + 2);
						}
						if (isOldFull(i + 1, j, N_x_old, N_y_old, old_full) && isOldFull(i, j - 1, N_x_old, N_y_old, old_full)) {
							mark(ii + 2, jj);
						}
						if (isOldFull(i - 1, j, N_x_old, N_y_old, old_full)) {
							mark(ii, jj + 1);
						}
						if (isOldFull(i + 1, j, N_x_old, N_y_old, old_full)) {
							mark(ii + 2, jj + 1);
						}
						if (isOldFull(i, j - 1, N_x_old, N_y_old, old_full)) {
							mark(ii + 1, jj);
						}
						if (isOldFull(i, j + 1, N_x_old, N_y_old, old_full)) {
							mark(ii + 1, jj + 2);
						}
					} else {
						keep(ii + 1, jj + 1, true);
						if (!isOldFull(i - 1, j, N_x_old, N_y_old, old_full) && !isOldFull(i, j - 1, N_x_old, N_y_old, old_full)) {
							keep(ii, jj, true);
						}
						if (!isOldFull(i - 1, j, N_x_old, N_y_old, old_full) && !isOldFull(i, j + 1, N_x_old, N_y_old, old_full)) {
							keep(ii, jj + 2, true);
						}
						if (!isOldFull(i + 1, j, N_x_old, N_y_old, old_full) && !isOldFull(i, j + 1, N_x_old, N_y_old, old_full)) {
							keep(ii + 2, jj + 2, true);
						}
						if (!isOldFull(i + 1, j, N_x_old, N_y_old, old_full) && !isOldFull(i, j - 1, N_x_old, N_y_old, old_full)) {
							keep(ii + 2, jj, true);
						}
						if (!isOldFull(i - 1, j, N_x_old, N_y_old, old_full)) {
							keep(ii, jj + 1, true);
						}
						if (!isOldFull(i + 1, j, N_x_old, N_y_old, old_full)) {
							keep(ii + 2, jj + 1, true);
						}
						if (!isOldFull(i, j - 1, N_x_old, N_y_old, old_full)) {
							keep(ii + 1, jj, true);
						}
						if (!isOldFull(i, j + 1, N_x_old, N_y_old, old_full)) {
							keep(ii + 1, jj + 2, true);
						}
					}
				}
			}

			if (hierarchyCntr <= hierarchy) {
				selectRandomCells();
			}

			hierarchyCntr += 1;
		}

		// output polygon

		// determine boundaries between marked and unmarked cells
		edges = makeEMatrix(nX + 1, nY + 1);

		if (diagonal) {
			top = makeBMatrix(nX, nY, false);
			bot = makeBMatrix(nX, nY, false);
			lft = makeBMatrix(nX, nY, false);
			rgt = makeBMatrix(nX, nY, false);
			for (i = 0; i < nX; ++i) {
				for (j = 0; j < nY; ++j) {
					if (isFull(i, j)) {
						setBot(i, j, true);
						setTop(i, j, true);
						setLft(i, j, true);
						setRgt(i, j, true);
					}
				}
			}
			for (i = 0; i <= nX; ++i) {
				for (j = 0; j <= nY; ++j) {
					if (isCompletelyFull(i, j)) {
						if (isFull(i - 1, j) && isFull(i - 1, j + 1) && isFull(i, j + 1) && isCompletelyEmpty(i, j - 1)
								&& isCompletelyEmpty(i + 1, j)) {
							setBot(i, j, false);
							setRgt(i, j, false);
						} else if (isFull(i + 1, j) && isFull(i + 1, j + 1) && isFull(i, j + 1) && isCompletelyEmpty(i, j - 1)
								&& isCompletelyEmpty(i - 1, j)) {
							setBot(i, j, false);
							setLft(i, j, false);
						} else if (isFull(i + 1, j) && isFull(i + 1, j - 1) && isFull(i, j - 1) && isCompletelyEmpty(i, j + 1)
								&& isCompletelyEmpty(i - 1, j)) {
							setTop(i, j, false);
							setLft(i, j, false);
						} else if (isFull(i - 1, j) && isFull(i - 1, j - 1) && isFull(i, j - 1) && isCompletelyEmpty(i, j + 1)
								&& isCompletelyEmpty(i + 1, j)) {
							setTop(i, j, false);
							setRgt(i, j, false);
						}
					}
				}
			}
			for (i = 0; i <= nX; ++i) {
				for (j = 0; j <= nY; ++j) {
					if (isCompletelyFull(i, j)) {
						if (isCompletelyEmpty(i, j - 1)) {
							storeEdge(i, j, i + 1, j);
						}
						if (isCompletelyEmpty(i - 1, j)) {
							storeEdge(i, j, i, j + 1);
						}
					} else if (isCompletelyEmpty(i, j)) {
						if (isCompletelyFull(i, j - 1)) {
							storeEdge(i, j, i + 1, j);
						}
						if (isCompletelyFull(i - 1, j)) {
							storeEdge(i, j + 1, i, j);

						}
					} else if (isTopFull(i, j)) {
						if (isRgtFull(i, j)) {
							storeEdge(i, j + 1, i + 1, j);
						} else if (isLftFull(i, j)) {
							storeEdge(i, j, i + 1, j + 1);
						}
					} else if (isBotFull(i, j)) {
						if (isLftFull(i, j)) {
							storeEdge(i, j + 1, i + 1, j);
						} else if (isRgtFull(i, j)) {
							storeEdge(i, j, i + 1, j + 1);
						}
					} else {
						System.err.format("cell[%s,%s]: neither top nor bottom is full!%n", i, j);
					}
				}
			}
		} else {
			for (i = 0; i <= nX; ++i) {
				for (j = 0; j <= nY; ++j) {
					if (isFull(i, j)) {
						if (!isFull(i, j - 1)) {
							storeEdge(i, j, i + 1, j);
						}
						if (!isFull(i - 1, j)) {
							storeEdge(i, j, i, j + 1);
						}
					} else {
						if (isFull(i, j - 1)) {
							storeEdge(i, j, i + 1, j);
						}
						if (isFull(i - 1, j)) {
							storeEdge(i, j, i, j + 1);
						}
					}
				}
			}
		}

		for (int i1 = 0; i1 <= nX; ++i1) {
			for (int j1 = 0; j1 <= nY; ++j1) {
				if (edges[i1][j1].i1 != NIL) {
					assert (edges[i1][j1].j1 != NIL);
					assert (edges[i1][j1].i2 != NIL);
					assert (edges[i1][j1].j2 != NIL);
					assert (((edges[edges[i1][j1].i1][edges[i1][j1].j1].i1 == i1) && (edges[edges[i1][j1].i1][edges[i1][j1].j1].j1 == j1))
							|| ((edges[edges[i1][j1].i1][edges[i1][j1].j1].i2 == i1)
									&& (edges[edges[i1][j1].i1][edges[i1][j1].j1].j2 == j1)));
				}
			}
		}

		final List<List<double[]>> rings = new ArrayList<>();
		int loopCntr = 0;
		for (int i1 = 0; i1 <= nX; ++i1) {
			for (int j1 = 0; j1 <= nY; ++j1) {
				if (edges[i1][j1].i1 != NIL) {
					assert (edges[i1][j1].j1 != NIL);
					assert (edges[i1][j1].i2 != NIL);
					assert (edges[i1][j1].j2 != NIL);
					assert (((edges[edges[i1][j1].i1][edges[i1][j1].j1].i1 == i1) && (edges[edges[i1][j1].i1][edges[i1][j1].j1].j1 == j1))
							|| ((edges[edges[i1][j1].i1][edges[i1][j1].j1].i2 == i1)
									&& (edges[edges[i1][j1].i1][edges[i1][j1].j1].j2 == j1)));
				}
				if (!getVis(i1, j1)) {
					rings.add(makePolygon(i1, j1, aligned, perturb, loopCntr, diagonal, smooth));
					++loopCntr;
				}
			}
		}

		return rings;
	}

	private void initializeGrid() {
		full = makeBMatrix(nX, nY, false);
		keep = makeBMatrix(nX, nY, false);
	}

	private List<double[]> makePolygon(int i1, int j1, final boolean aligned, final boolean perturb, final int loop_cntr,
			final boolean diagonal, int smooth) {
		int number = 0, i0, j0;
		int i2, j2, sum = 0;
		VertexNode zero = new VertexNode(0, 0);
		final Coord p = new Coord(0, 0);
		List<Coord> oldPnts;

		assert ((i1 >= 0) && (i1 <= nX) && (j1 >= 0) && (j1 <= nY));

		// count the number of edges of this loop
		i0 = i1;
		j0 = j1;
		numVertices = 0;

		do {
			storeVertex(i1, j1);
			setVis(i1, j1, true);
			i2 = getStartI(i1, j1);
			j2 = getStartJ(i1, j1);
			assert ((i2 >= 0) && (i2 <= nX) && (j2 >= 0) && (j2 <= nY));
			if (getVis(i2, j2)) {
				i2 = getEndI(i1, j1);
				j2 = getEndJ(i1, j1);
				if (getVis(i2, j2)) {
					i2 = i0;
					j2 = j0;
				}
				assert ((i2 >= 0) && (i2 <= nX) && (j2 >= 0) && (j2 <= nY));
			}
			i1 = i2;
			j1 = j2;
		} while (!getVis(i1, j1));
		assert ((i0 == i1) && (j0 == j1));
		storeVertex(i0, j0);
		storeVertex(i0, j0);
		number = numVertices - 1;

		int i = 0;
		int j = 1;
		int k = 2;
		while (k < number) {
			while ((k < number) && (((vertices[i].i1 == vertices[j].i1) && (vertices[i].i1 == vertices[k].i1))
					|| ((vertices[i].j1 == vertices[j].j1) && (vertices[i].j1 == vertices[k].j1)))) {
				j = k;
				++k;
			}
			++i;
			vertices[i].i1 = vertices[j].i1;
			vertices[i].j1 = vertices[j].j1;
			++j;
			++k;
		}
		++i;
		vertices[i].i1 = vertices[j].i1;
		vertices[i].j1 = vertices[j].j1;

		if (diagonal) {
			number = i + 1;
			i = 0;
			j = 1;
			k = 2;
			while (k < number) {
				while ((k < number)
						&& ((vertices[i].i1 != vertices[j].i1) && (vertices[i].j1 != vertices[j].j1) && (vertices[i].i1 != vertices[k].i1)
								&& (vertices[i].j1 != vertices[k].j1))
						&& (((vertices[i].i1 - vertices[j].i1) * (vertices[i].j1 - vertices[k].j1)) == ((vertices[i].j1 - vertices[j].j1)
								* (vertices[i].i1 - vertices[k].i1)))) {
					j = k;
					++k;
				}
				++i;
				vertices[i].i1 = vertices[j].i1;
				vertices[i].j1 = vertices[j].j1;
				++j;
				++k;
			}
			++i;
			vertices[i].i1 = vertices[j].i1;
			vertices[i].j1 = vertices[j].j1;
		}

		if ((vertices[i - 1].i1 == vertices[i].i1) && (vertices[i - 1].j1 == vertices[i].j1)) {
			number = i;
		} else if ((vertices[0].i1 == vertices[i].i1) && (vertices[0].j1 == vertices[i].j1)) {
			number = i + 1;
		} else {
			number = i;
		}

		for (i = 1; i < number; ++i) {
			sum += det2D(vertices[i - 1], vertices[i], zero);
		}

		if (((loop_cntr == 0) && (sum < 0)) || ((loop_cntr > 0) && (sum > 0))) {
			i = 1;
			j = number - 2;
			while (i < j) {
				zero = vertices[i]; // swap
				vertices[i] = vertices[j];
				vertices[j] = zero;
				++i;
				--j;
			}
		}

		final List<double[]> ring = new ArrayList<>(number);

		if (aligned) {
			for (int l = 0; l < number; l++) {
				ring.add(new double[] { vertices[l].i1, vertices[l].j1 });
			}
		} else {
			pnts = new ArrayList<>((number - 1) * (smooth + 1) + 2);
			if (perturb) {
				--number;
				for (i = 0; i < number; ++i) {
					p.x = vertices[i].i1 + perturbation();
					p.y = vertices[i].j1 + perturbation();
					storePnt(p);
				}
			} else {
				p.x = vertices[0].i1; // could just instantiate p
				p.y = vertices[0].j1;
				storePnt(p);
				number -= 2;
				for (i = 1; i < number; ++i) {
					if (vertices[i].i1 == vertices[i - 1].i1) {
						p.y = vertices[i].j1 + perturbation();
					} else {
						p.x = vertices[i].i1 + perturbation();
					}
					storePnt(p);
				}
				i = number;
				if (vertices[i].i1 == vertices[i - 1].i1) {
					p.y = vertices[i].j1;
				} else {
					p.x = vertices[i].i1;
				}
				storePnt(p);
			}
			storePnt(pnts.get(0)); // close ring (unperturbed)

			while (smooth > 0) {
				oldPnts = pnts;
				pnts = new ArrayList<>(oldPnts.size() * 2);
				for (i = 1; i < oldPnts.size(); ++i) {
					p.x = (3.0 * oldPnts.get(i - 1).x + oldPnts.get(i).x) / 4.0;
					p.y = (3.0 * oldPnts.get(i - 1).y + oldPnts.get(i).y) / 4.0;
					storePnt(p);
					p.x = (oldPnts.get(i - 1).x + 3.0 * oldPnts.get(i).x) / 4.0;
					p.y = (oldPnts.get(i - 1).y + 3.0 * oldPnts.get(i).y) / 4.0;
					storePnt(p);
				}
				p.x = (3.0 * oldPnts.get(0).x + oldPnts.get(1).x) / 4.0;
				p.y = (3.0 * oldPnts.get(0).y + oldPnts.get(1).y) / 4.0;
				storePnt(p);
				--smooth;
			}

			for (int l = 0; l < pnts.size(); l++) {
				ring.add(new double[] { pnts.get(l).x, pnts.get(l).y });
			}
		}
		pnts = null;

		return ring;
	}

	private int uniformRandom(final int m) {
		return rand.nextInt(m);
	}

	private int convert(final int i, final int j) {
		return i * nY + j;
	}

	private static double det2D(final VertexNode u, final VertexNode v, final VertexNode w) {
		return (((u).i1 - (v).i1) * ((v).j1 - (w).j1) + ((v).j1 - (u).j1) * ((v).i1 - (w).i1));
	}

	private void setTop(final int I, final int J, final boolean W) {
		assert (I >= 0 && I < nX && J >= 0 && J < nY);
		top[I][J] = W;
	}

	private void setBot(final int I, final int J, final boolean W) {
		assert (I >= 0 && I < nX && J >= 0 && J < nY);
		bot[I][J] = W;
	}

	private void setLft(final int I, final int J, final boolean W) {
		assert (I >= 0 && I < nX && J >= 0 && J < nY);
		lft[I][J] = W;
	}

	private void setRgt(final int I, final int J, final boolean W) {
		assert (I >= 0 && I < nX && J >= 0 && J < nY);
		rgt[I][J] = W;
	}

	private boolean isTopFull(final int I, final int J) {
		return top[I][J];
	}

	private boolean isBotFull(final int i, final int j) {
		return bot[i][j];
	}

	private boolean isLftFull(final int i, final int j) {
		return lft[i][j];
	}

	private boolean isRgtFull(final int i, final int j) {
		return rgt[i][j];
	}

	private int getStartI(final int i, final int j) {
		assert (i >= 0 && i <= nX && j >= 0 && j <= nY);
		return edges[i][j].i1;
	}

	private int getStartJ(final int i, final int j) {
		assert (i >= 0 && i <= nX && j >= 0 && j <= nY);
		return edges[i][j].j1;
	}

	private int getEndI(final int i, final int j) {
		assert (i >= 0 && i <= nX && j >= 0 && j <= nY);
		return edges[i][j].i2;
	}

	private int getEndJ(final int i, final int j) {
		assert (i >= 0 && i <= nX && j >= 0 && j <= nY);
		return edges[i][j].j2;
	}

	private boolean getVis(final int i, final int j) {
		assert (i >= 0 && i <= nX && j >= 0 && j <= nY);
		return edges[i][j].vis;
	}

	private void setVis(final int i, final int j, final boolean b) {
		assert (i >= 0 && i <= nX && j >= 0 && j <= nY);
		edges[i][j].vis = b;
	}

	private boolean[][] makeBMatrix(final int Nx, final int Ny, final boolean value) {
		final boolean[][] matrix = new boolean[Nx][Ny];

		for (int i = 0; i < Nx; ++i) {
			for (int j = 0; j < Ny; ++j) {
				matrix[i][j] = value;
			}
		}

		return matrix;
	}

	private EdgeNode[][] makeEMatrix(final int Nx, final int Ny) {
		int i, j;
		final EdgeNode[][] matrix = new EdgeNode[Nx][];

		for (i = 0; i < Nx; ++i) {
			matrix[i] = new EdgeNode[Ny];
			for (j = 0; j < Ny; ++j) {
				matrix[i][j] = new EdgeNode(NIL, NIL, NIL, NIL, true);
			}
		}

		return matrix;
	}

	private double perturbation() {
		final double max = 400 / 899d;
		return rand.nextDouble(-max, max);
	}

	private void storePnt(final Coord P) {
		pnts.add(new Coord(P.x, P.y));
	}

	private void storeEdge(final int i1, final int j1, final int i2, final int j2) {
		if ((i1 >= 0) && (j1 >= 0) && (i1 <= nX) && (j1 <= nY) && (i2 >= 0) && (j2 >= 0) && (i2 <= nX) && (j2 <= nY)) {
			setVis(i1, j1, false);
			if (edges[i1][j1].i1 == NIL) {
				assert edges[i1][j1].j1 == NIL;
				edges[i1][j1].i1 = i2;
				edges[i1][j1].j1 = j2;
			} else {
				assert edges[i1][j1].j1 != NIL;
				assert edges[i1][j1].i2 == NIL;
				assert edges[i1][j1].j2 == NIL;
				edges[i1][j1].i2 = i2;
				edges[i1][j1].j2 = j2;
			}
			setVis(i2, j2, false);
			if (edges[i2][j2].i1 == NIL) {
				assert edges[i2][j2].j1 == NIL;
				edges[i2][j2].i1 = i1;
				edges[i2][j2].j1 = j1;
			} else {
				assert edges[i2][j2].j1 != NIL;
				assert edges[i2][j2].i2 == NIL;
				assert edges[i2][j2].j2 == NIL;
				edges[i2][j2].i2 = i1;
				edges[i2][j2].j2 = j1;
			}
		} else {
			System.err.println("StoreEdge(): index out of bounds!");
			System.err.printf("(%d,%d) <--> (%d,%d)%n", i1, j1, i2, j2);
		}
	}

	private void storeVertex(final int i1, final int j1) {
		if (numVertices >= maxNumVertices) {
			maxNumVertices += 500;
			vertices = Arrays.copyOf(vertices, maxNumVertices);
		}
		vertices[numVertices] = new VertexNode(i1, j1);
		numVertices++;
	}

	private boolean isFull(final int i, final int j) {
		if (i >= 0 && i < nX && j >= 0 && j < nY) {
			return full[i][j];
		} else {
			return false;
		}
	}

	private boolean isOldFull(final int i, final int j, final int N_x_old, final int N_y_old, final boolean[][] old_full) {
		if (i >= 0 && i < N_x_old && j >= 0 && j < N_y_old) {
			return old_full[i][j];
		} else {
			return false;
		}
	}

	private boolean isCompletelyFull(final int i, final int j) {
		if (i >= 0 && i < nX && j >= 0 && j < nY) {
			return (top[i][j] && bot[i][j] && lft[i][j] && rgt[i][j]);
		} else {
			return false;
		}
	}

	private boolean isCompletelyEmpty(final int i, final int j) {
		if (i >= 0 && i < nX && j >= 0 && j < nY) {
			return (!(top[i][j] || bot[i][j] || lft[i][j] || rgt[i][j]));
		} else {
			return true;
		}
	}

	private boolean toBeKept(final int i, final int j) {
		if ((i >= 0) && (i < nX) && (j >= 0) && (j < nY)) {
			return keep[i][j];
		} else {
			return true;
		}
	}

	private boolean isPossible(final int i, final int j) {
		int c = 0;

		if ((i < 0) || (i >= nX) || (j < 0) || (j >= nY)) {
			return false;
		}

		if (toBeKept(i, j) || isFull(i, j)) {
			return false;
		}

		if (isFull(i - 1, j + 1)) {
			if (!(isFull(i - 1, j) || isFull(i, j + 1))) {
				return false;
			}
		}
		if (isFull(i - 1, j - 1)) {
			if (!(isFull(i - 1, j) || isFull(i, j - 1))) {
				return false;
			}
		}
		if (isFull(i + 1, j - 1)) {
			if (!(isFull(i + 1, j) || isFull(i, j - 1))) {
				return false;
			}
		}
		if (isFull(i + 1, j + 1)) {
			if (!(isFull(i + 1, j) || isFull(i, j + 1))) {
				return false;
			}
		}

		if (holes) {
			c = uniformRandom(30);
		}

		if (c != 0) {
			if (isFull(i, j - 1) && isFull(i, j + 1)) {
				if (!(isFull(i + 1, j) || isFull(i - 1, j))) {
					return false;
				}
			}
			if (isFull(i - 1, j) && isFull(i + 1, j)) {
				if (!(isFull(i, j - 1) || isFull(i, j + 1))) {
					return false;
				}
			}
		}

		return true;
	}

	private void storeCandidate(final int i, final int j) {
		int k;

		if ((i >= 0) && (i < nX) && (j >= 0) && (j < nY)) {
			if (!(toBeKept(i, j) && isFull(i, j))) {
				k = convert(i, j);
				if (numCandidates >= candidates.size()) {
					candidates.add(k);
				} else {
					candidates.set(numCandidates, k);
				}
				numCandidates++;
			}
		}
	}

	private void storeKeepCandidate(final int i, final int j) {
		int k;

		if ((i >= 0) && (i < nX) && (j >= 0) && (j < nY)) {
			if (!(toBeKept(i, j) && isFull(i, j))) {
				k = convert(i, j);
				if (numKeepCandidates >= keepCandidates.size()) {
					keepCandidates.add(k);
				} else {
					keepCandidates.set(numKeepCandidates, k);
				}
				numKeepCandidates++;
			}
		}
	}

	private void mark(final int i, final int j) {
		assert ((i >= 0) && (i < nX) && (j >= 0) && (j < nY));
		full[i][j] = true;
		numCells++;

		if (i > 0) {
			storeCandidate(i - 1, j);
		}
		if (j > 0) {
			storeCandidate(i, j - 1);
		}
		if (i < (nX - 1)) {
			storeCandidate(i + 1, j);
		}
		if (j < (nY - 1)) {
			storeCandidate(i, j + 1);
		}
	}

	private void keep(final int i, final int j, final boolean storeCandidates) {
		if ((i >= 0) && (j >= 0) && (i < nX) && (j < nY)) {
			keep[i][j] = true;
			numKeep++;

			if (storeCandidates) {
				if (i > 0) {
					storeKeepCandidate(i - 1, j);
				}
				if (j > 0) {
					storeKeepCandidate(i, j - 1);
				}
				if (i < (nX - 1)) {
					storeKeepCandidate(i + 1, j);
				}
				if (j < (nY - 1)) {
					storeKeepCandidate(i, j + 1);
				}
			}
		}
	}

	private void selectRandomCells() {
		int k, m, mm, i = 0, j, c = 0;

		while ((numCells < maxCells) && (numCandidates > 0)) {
			c = uniformRandom(numCandidates);
			k = candidates.get(c);
			--numCandidates;
			candidates.set(c, candidates.get(numCandidates));
			i = k / nY;
			j = k - i * nY;
			if (isPossible(i, j)) {
				mark(i, j);
			}
			if ((numKeepCandidates > 0) && (numKeep < maxKeep)) {
				if ((maxCells - numCells - 1) > 0) {
					mm = 1 + 2 * (maxKeep - numKeep) / (maxCells - numCells - 1);
				} else {
					mm = 1 + 2 * (maxKeep - numKeep);
				}
				m = 0;
				while ((m < mm) && (numKeepCandidates > 0)) {
					c = uniformRandom(numKeepCandidates);
					k = keepCandidates.get(c);
					--numKeepCandidates;
					keepCandidates.set(c, keepCandidates.get(numKeepCandidates));
					i = k / nY;
					j = k - i * nY;
					if (!isFull(i, j)) {
						keep(i, j, true);
						++m;
					}
				}
			}
		}
	}

	private static class Coord {
		double x;
		double y;

		private Coord(final double x, final double y) {
			this.x = x;
			this.y = y;
		}

	}

	private static class EdgeNode {
		int i1;
		int j1;
		int i2;
		int j2;
		boolean vis;

		private EdgeNode(final int i1, final int j1, final int i2, final int j2, final boolean vis) {
			this.i1 = i1;
			this.j1 = j1;
			this.i2 = i2;
			this.j2 = j2;
			this.vis = vis;
		}

	}

	private static class VertexNode {
		int i1;
		int j1;

		private VertexNode(final int i1, final int j1) {
			this.i1 = i1;
			this.j1 = j1;
		}

	}

}

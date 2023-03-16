package micycle.srpg;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * SRPolygonGenerator - Super Random Polygon Generator
 * <p>
 * SRPolygonGenerator generates simply-connected and multiply-connected polygons by means of a
 * regular grid that consists of square cells. Given two integer values, a and
 * b, SRPolygonGenerator generates a grid of size a times b.
 * <P>
 * By default SRPolygonGenerator then generates orthogonal polygons on this grid. An
 * additional parameter p, between zero and one, leads to a smaller or larger
 * number of vertices in the produced polygon. SRPolygonGenerator is able to produce octagonal
 * polygons by cutting off corners with ±45° diagonals during the construction.
 * Cutting corners repeatedly, without the diagonal restriction, yields an
 * approximation of a smooth free-form curve. Additionally, SRPolygonGenerator can apply
 * perturbations in order to generate polygons with axes-parallel edges whose
 * vertices do not lie on a grid, or to generate polygons whose edges (in
 * general) are not parallel to the coordinate axes.
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

	// https://github.com/cgalab/genpoly-srpg/blob/master/srpg.c
	// https://sbgdb.cs.sbg.ac.at/classes/polygons/
	// https://www1.pub.informatik.uni-wuerzburg.de/eurocg2020/data/uploads/papers/eurocg20_paper_75.pdf

	private static final int NIL = -1;

	private static class Coord {
		double x;
		double y;

		public Coord(double x, double y) {
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

		public EdgeNode(int i1, int j1, int i2, int j2, boolean vis) {
			this.i1 = i1;
			this.j1 = j1;
			this.i2 = i2;
			this.j2 = j2;
			this.vis = vis;
		}

		public EdgeNode() {
		}

	}

	private static class VertexNode {
		int i1;
		int j1;

		public VertexNode(int i1, int j1) {
			this.i1 = i1;
			this.j1 = j1;
		}

	}

	private boolean[][] full = null;
	private boolean[][] keep = null;
	private int N = 0;

	private boolean[][] top = null;
	private boolean[][] bot = null;
	private boolean[][] lft = null;
	private boolean[][] rgt = null;

	private int num_candidates = 0;
	private int max_num_candidates = 0;
	private int num_keep_candidates = 0;
	private int max_num_keep_candidates = 0;

	private EdgeNode[][] edges = null;

	private int num_vertices = 0;
	private int max_num_vertices = 0;

	int state = 1337;

	int xorshift() {
		int x = state;
		x ^= x << 13;
		x ^= x >> 17;
		x ^= x << 5;
		state = Math.abs(x);
		return state;
	}

	private int uniformRandom(int m) { // NOTE C ARG
		if (m == 0) {
			m++;
		}
		return xorshift() % m;
//		return rand.nextInt(m);
	}

	private int convert(int i, int j) {
		return i * N_y + j;
	}

	private int invert(int k, int i) {
		i = k / N_y; // TODO has side effects / should mutate?
		return k - i * N_y;
	}

	private double det2D(VertexNode u, VertexNode v, VertexNode w) {
		return (((u).i1 - (v).i1) * ((v).j1 - (w).j1) + ((v).j1 - (u).j1) * ((v).i1 - (w).i1));
	}

//	@Deprecated
//	private void swap(VertexNode i1, VertexNode i2, VertexNode i) { // TODO does not work in java -- replace?!
//		i = i1;
//		i1 = i2;
//		i2 = i;
//	}

	private void swap(VertexNode[] vertices, int i, int j) {
		// TODO replace inline swap with this (once confirmed working)
		VertexNode temp = vertices[i];
		vertices[i] = vertices[j];
		vertices[j] = temp;
	}

	private void setTop(int I, int J, boolean W) {
		assert (I >= 0 && I < N_x && J >= 0 && J < N_y);
		top[I][J] = W;
	}

	private void setBot(int I, int J, boolean W) {
		assert (I >= 0 && I < N_x && J >= 0 && J < N_y);
		bot[I][J] = W;
	}

	private void setLft(int I, int J, boolean W) {
		assert (I >= 0 && I < N_x && J >= 0 && J < N_y);
		lft[I][J] = W;
	}

	private void setRgt(int I, int J, boolean W) {
		assert (I >= 0 && I < N_x && J >= 0 && J < N_y);
		rgt[I][J] = W;
	}

	private boolean isTopFull(int I, int J) {
		return top[I][J];
	}

	private boolean IsBotFull(int i, int j) {
		return bot[i][j];
	}

	private boolean IsLftFull(int i, int j) {
		return lft[i][j];
	}

	private boolean IsRgtFull(int i, int j) {
		return rgt[i][j];
	}

	private int GetStartI(int i, int j) {
		assert (i >= 0 && i <= N_x && j >= 0 && j <= N_y);
		return edges[i][j].i1;
	}

	private int GetStartJ(int i, int j) {
		assert (i >= 0 && i <= N_x && j >= 0 && j <= N_y);
		return edges[i][j].j1;
	}

	private int GetEndI(int i, int j) {
		assert (i >= 0 && i <= N_x && j >= 0 && j <= N_y);
		return edges[i][j].i2;
	}

	private int GetEndJ(int i, int j) {
		assert (i >= 0 && i <= N_x && j >= 0 && j <= N_y);
		return edges[i][j].j2;
	}

	private boolean GetVis(int i, int j) {
		assert (i >= 0 && i <= N_x && j >= 0 && j <= N_y);
		return edges[i][j].vis;
	}

	private void SetVis(int i, int j, boolean b) {
		assert (i >= 0 && i <= N_x && j >= 0 && j <= N_y);
		edges[i][j].vis = b;
	}

	public boolean[][] makeBMatrix(int Nx, int Ny, boolean value) {
		boolean[][] matrix = new boolean[Nx][Ny];

		for (int i = 0; i < Nx; ++i) {
			for (int j = 0; j < Ny; ++j) {
				matrix[i][j] = value;
			}
		}

		return matrix;
	}

	public EdgeNode[][] MakeEMatrix(int Nx, int Ny) {
		int i, j;
		EdgeNode[][] matrix = new EdgeNode[Nx][];

		for (i = 0; i < Nx; ++i) {
			matrix[i] = new EdgeNode[Ny];
			for (j = 0; j < Ny; ++j) {
				matrix[i][j] = new EdgeNode(NIL, NIL, NIL, NIL, true);
//				matrix[i][j].i1 = NIL;
//				matrix[i][j].j1 = NIL;
//				matrix[i][j].i2 = NIL;
//				matrix[i][j].j2 = NIL;
//				matrix[i][j].vis = true;
			}
		}

		return matrix;
	}

	double Perturbation() {
		int c = 0;

		c = uniformRandom(800001);
		c -= 400000;

		return ((c) / 899000.0);
	}

	private void storePnt(Coord P) {
		if (num_pnts >= max_num_pnts) {
			max_num_pnts += 1001;
			pnts = Arrays.copyOf(pnts, max_num_pnts);
		}

		pnts[num_pnts] = new Coord(P.x, P.y);
		num_pnts++;
	}

	void storeEdge(int i1, int j1, int i2, int j2) {
		if ((i1 >= 0) && (j1 >= 0) && (i1 <= N_x) && (j1 <= N_y) && (i2 >= 0) && (j2 >= 0) && (i2 <= N_x) && (j2 <= N_y)) {
			SetVis(i1, j1, false);
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
			SetVis(i2, j2, false);
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
			System.err.printf("             (%d,%d) <--> (%d,%d)\n", i1, j1, i2, j2);
		}
	}

	public void storeVertex(int i1, int j1) {
		if (num_vertices >= max_num_vertices) {
			max_num_vertices += 1001;
			vertices = Arrays.copyOf(vertices, max_num_vertices);
		}
		vertices[num_vertices] = new VertexNode(i1, j1);
		num_vertices++;
	}

	public boolean IsFull(int i, int j) {
		if (i >= 0 && i < N_x && j >= 0 && j < N_y) {
			return full[i][j];
		} else {
			return false;
		}
	}

	public boolean IsOldFull(int i, int j, int N_x_old, int N_y_old, boolean[][] old_full) {
		if (i >= 0 && i < N_x_old && j >= 0 && j < N_y_old) {
			return old_full[i][j];
		} else {
			return false;
		}
	}

	public boolean isCompletelyFull(int i, int j) {
		if (i >= 0 && i < N_x && j >= 0 && j < N_y) {
			return (top[i][j] && bot[i][j] && lft[i][j] && rgt[i][j]);
		} else {
			return false;
		}
	}

	public boolean isCompletelyEmpty(int i, int j) {
		if (i >= 0 && i < N_x && j >= 0 && j < N_y) {
			return (!(top[i][j] || bot[i][j] || lft[i][j] || rgt[i][j]));
		} else {
			return true;
		}
	}

	boolean ToBeKept(int i, int j) {
		if ((i >= 0) && (i < N_x) && (j >= 0) && (j < N_y)) {
			return keep[i][j];
		} else {
			return true;
		}
	}

	boolean IsPossible(int i, int j) {
		int c = 0;

		if ((i < 0) || (i >= N_x) || (j < 0) || (j >= N_y)) {
			return false;
		}

		if (ToBeKept(i, j) || IsFull(i, j)) {
			return false;
		}

		if (IsFull(i - 1, j + 1)) {
			if (!(IsFull(i - 1, j) || IsFull(i, j + 1))) {
				return false;
			}
		}
		if (IsFull(i - 1, j - 1)) {
			if (!(IsFull(i - 1, j) || IsFull(i, j - 1))) {
				return false;
			}
		}
		if (IsFull(i + 1, j - 1)) {
			if (!(IsFull(i + 1, j) || IsFull(i, j - 1))) {
				return false;
			}
		}
		if (IsFull(i + 1, j + 1)) {
			if (!(IsFull(i + 1, j) || IsFull(i, j + 1))) {
				return false;
			}
		}

		if (holes) {
			c = uniformRandom(30);
		}

		if (c != 0) {
			if (IsFull(i, j - 1) && IsFull(i, j + 1)) {
				if (!(IsFull(i + 1, j) || IsFull(i - 1, j))) {
					return false;
				}
			}
			if (IsFull(i - 1, j) && IsFull(i + 1, j)) {
				if (!(IsFull(i, j - 1) || IsFull(i, j + 1))) {
					return false;
				}
			}
		}

		return true;
	}

	void StoreCandidate(int i, int j) {
		int k;

		if ((i >= 0) && (i < N_x) && (j >= 0) && (j < N_y)) {
			if (!(ToBeKept(i, j) && IsFull(i, j))) {
				k = convert(i, j); // NOTE arg change
				if (num_candidates >= max_num_candidates) {
					max_num_candidates += 1001;
					candidates = Arrays.copyOf(candidates, max_num_candidates);
				}
				candidates[num_candidates] = k;
				++num_candidates;
			}
		}
	}

	void StoreKeepCandidate(int i, int j) {
		int k;

		if ((i >= 0) && (i < N_x) && (j >= 0) && (j < N_y)) {
			if (!(ToBeKept(i, j) && IsFull(i, j))) {
				k = convert(i, j);
				if (num_keep_candidates >= max_num_keep_candidates) {
					max_num_keep_candidates += 1001;
					keep_candidates = Arrays.copyOf(keep_candidates, max_num_keep_candidates);
				}
				keep_candidates[num_keep_candidates] = k;
				++num_keep_candidates;
			}
		}
	}

	void Mark(int i, int j, int[] num_cells) { // NOTE C STYLE INPUT int*
		assert ((i >= 0) && (i < N_x) && (j >= 0) && (j < N_y));
		full[i][j] = true;
		num_cells[0]++;

		if (i > 0) {
			StoreCandidate(i - 1, j);
		}
		if (j > 0) {
			StoreCandidate(i, j - 1);
		}
		if (i < (N_x - 1)) {
			StoreCandidate(i + 1, j);
		}
		if (j < (N_y - 1)) {
			StoreCandidate(i, j + 1);
		}
	}

	void Keep(int i, int j, boolean store_candidates, int[] num_keep) { // NOTE C STYLE INPUT int*
		if ((i >= 0) && (j >= 0) && (i < N_x) && (j < N_y)) {
			keep[i][j] = true;
			num_keep[0]++;

			if (store_candidates) {
				if (i > 0) {
					StoreKeepCandidate(i - 1, j);
				}
				if (j > 0) {
					StoreKeepCandidate(i, j - 1);
				}
				if (i < (N_x - 1)) {
					StoreKeepCandidate(i + 1, j);
				}
				if (j < (N_y - 1)) {
					StoreKeepCandidate(i, j + 1);
				}
			}
		}
	}

	int num_pnts = 0, max_num_pnts; // NOTE these are not defined outside method in C version

	List<double[]> OutputLoop(PrintStream output, int i1, int j1, int[] total_number, boolean aligned, boolean perturb, int loop_cntr,
			boolean diagonal, int smooth) {
		int number = 0, i0, j0;
//int e_comp(const void *e1, const void *e2); // NOTE WTF?
		int i, j, k, i2, j2, sum = 0;
		VertexNode zero = new VertexNode(0, 0);
		Coord p = new Coord(0, 0); // NOTE CHECK
		pnts = null; // NOTE CHECK
		Coord[] old_pnts = null; // NOTE CHECK
		num_pnts = 0;
		max_num_pnts = 0;
		int old_num_pnts;

		assert ((i1 >= 0) && (i1 <= N_x) && (j1 >= 0) && (j1 <= N_y));

		/*                                                                        */
		/* count the number of edges of this loop */
		/*                                                                        */
		i0 = i1;
		j0 = j1;
		num_vertices = 0;

		do {
			storeVertex(i1, j1);
			SetVis(i1, j1, true);
			i2 = GetStartI(i1, j1);
			j2 = GetStartJ(i1, j1);
			assert ((i2 >= 0) && (i2 <= N_x) && (j2 >= 0) && (j2 <= N_y));
			if (GetVis(i2, j2)) {
				i2 = GetEndI(i1, j1);
				j2 = GetEndJ(i1, j1);
				if (GetVis(i2, j2)) {
					i2 = i0;
					j2 = j0;
				}
				assert ((i2 >= 0) && (i2 <= N_x) && (j2 >= 0) && (j2 <= N_y));
			}
			i1 = i2;
			j1 = j2;
		} while (!GetVis(i1, j1));
		assert ((i0 == i1) && (j0 == j1));
		storeVertex(i0, j0);
		storeVertex(i0, j0);
		number = num_vertices - 1;

		i = 0;
		j = 1;
		k = 2;
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
				// NOTE inlined swap() method here
//				swap(vertices[i], vertices[j], zero);
				zero = vertices[i];
				vertices[i] = vertices[j];
				vertices[j] = zero;
				++i;
				--j;
			}
		}

		List<double[]> ring = new ArrayList<>(vertices.length);

		if (aligned) {
			total_number[0] += number; // *total_number += number;
//			output.printf("%d%n", number);
//			output.printf("%d %d\n", vertices[0].i1, vertices[0].j1);
//			for (i = 1; i < number; ++i) {
//				output.printf("%d %d\n", vertices[i].i1, vertices[i].j1);
//			}
//			output.printf("\n");

			for (int l = 0; l < vertices.length; l++) {
				ring.add(new double[] {vertices[l].i1, vertices[l].j1});
			}
		} else {
			max_num_pnts = (number - 1) * (smooth + 1) + 2;
			num_pnts = 0;
			pnts = new Coord[max_num_pnts]; // NOTE
			if (perturb) {
				--number;
				for (i = 0; i < number; ++i) {
					p.x = vertices[i].i1 + Perturbation();
					p.y = vertices[i].j1 + Perturbation();
					storePnt(p);
				}
				storePnt(pnts[0]); // close ring (unperturbed)
			} else {
				p.x = vertices[0].i1; // could just instantiate p
				p.y = vertices[0].j1;
				storePnt(p);
				number -= 2;
				for (i = 1; i < number; ++i) {
					if (vertices[i].i1 == vertices[i - 1].i1) {
						p.y = vertices[i].j1 + Perturbation();
					} else {
						p.x = vertices[i].i1 + Perturbation();
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
				storePnt(pnts[0]); // close ring (unperturbed)
			}

			while (smooth > 0) {
				old_pnts = pnts;
				old_num_pnts = num_pnts;
				max_num_pnts *= 2;
				num_pnts = 0;
				pnts = new Coord[max_num_pnts]; // NOTE
				for (i = 1; i < old_num_pnts; ++i) {
					p.x = (3.0 * old_pnts[i - 1].x + old_pnts[i].x) / 4.0;
					p.y = (3.0 * old_pnts[i - 1].y + old_pnts[i].y) / 4.0;
					storePnt(p);
					p.x = (old_pnts[i - 1].x + 3.0 * old_pnts[i].x) / 4.0;
					p.y = (old_pnts[i - 1].y + 3.0 * old_pnts[i].y) / 4.0;
					storePnt(p);
				}
				p.x = (3.0 * old_pnts[0].x + old_pnts[1].x) / 4.0;
				p.y = (3.0 * old_pnts[0].y + old_pnts[1].y) / 4.0;
				storePnt(p);
				--smooth;
			}

//			output.printf("%d\n", num_pnts);
//			output.printf("%s, %s\n", pnts[0].x, pnts[0].y);
//			for (i = 1; i < num_pnts; ++i) {
//				output.printf("%s %s\n", pnts[i].x, pnts[i].y);
//			}
//			output.printf("\n");

			for (int l = 0; l < number; l++) {
				ring.add(new double[] { pnts[l].x, pnts[l].y });
			}

			total_number[0] += num_pnts; // *total_number += num_pnts;
		}
		pnts = null;

		ring.add(ring.get(0));
		return ring;
	}

	void initializeGrid() {
		full = makeBMatrix(N_x, N_y, false);
		keep = makeBMatrix(N_x, N_y, false);
	}

	void SelectRandomCells(int[] num_cells, int max_cells, int[] num_keep, int max_keep) { // NOTE C-style arg passing int*
		int k, m, mm, i = 0, j, c = 0;

		while ((num_cells[0] < max_cells) && (num_candidates > 0)) {
			c = uniformRandom(num_candidates);
			k = candidates[c];
			--num_candidates;
			candidates[c] = candidates[num_candidates];
//			j = invert(k, i);
			i = k / N_y; // NOTE inlined invert()
			j = k - i * N_y; // NOTE inlined invert()
			if (IsPossible(i, j)) {
				Mark(i, j, num_cells);
			}
			if ((num_keep_candidates > 0) && (num_keep[0] < max_keep)) {
				if ((max_cells - num_cells[0] - 1) > 0) {
					mm = 1 + 2 * (max_keep - num_keep[0]) / (max_cells - num_cells[0] - 1);
				} else {
					mm = 1 + 2 * (max_keep - num_keep[0]);
				}
				m = 0;
				while ((m < mm) && (num_keep_candidates > 0)) {
					c = uniformRandom(num_keep_candidates);
					k = keep_candidates[c];
					--num_keep_candidates;
					keep_candidates[c] = keep_candidates[num_keep_candidates];
//					j = invert(k, i);
					i = k / N_y; // NOTE INLINED INVERT
					j = k - i * N_y; // NOTE INLINED INVERT
					if (!IsFull(i, j)) {
						Keep(i, j, true, num_keep);
						++m;
					}
				}
			}
		}
	}

	PrintStream output = System.out;
	int i = 0, j = 0, k, m, n, mm, nn;
	int max_keep, max_cells;
	int[] num_keep = new int[1], num_cells = new int[1]; // NOTE replicates c++ int* (pass by value)
	int[] total_number = new int[1];
	int hierarchy_cntr = 1;
	int loop_cntr, i1, j1;
	double keep_percent;
	int N_x_old, N_y_old, ii, jj;
	boolean[][] old_keep, old_full;

	// generates the polygon
	public List<List<double[]>> Compute(double mark_percent, boolean perturb, boolean aligned, int hierarchy, boolean diagonal, int smooth,
			String file_name) {
		// 1. allocate the grid
		N = N_x * N_y;
		initializeGrid();

		// 2. select cells to be kept
		keep_percent = (1.0 - mark_percent) * 0.95;
		max_keep = (int) (N * keep_percent);
		max_cells = (int) (N * mark_percent);

		m = N_y / 10;
		for (i = 0; i < N_x; ++i) {
			k = uniformRandom(m);
			for (j = 0; j <= k; ++j) {
				Keep(i, j, false, num_keep);
			}
			k = uniformRandom(m);
			for (j = N_y - 1; j >= N_y - k - 1; --j) {
				Keep(i, j, false, num_keep);
			}
		}
		for (j = 0; j < N_y; ++j) {
			k = uniformRandom(m);
			for (i = 0; i <= k; ++i) {
				Keep(i, j, false, num_keep);
			}
			k = uniformRandom(m);
			for (i = N_x - 1; i >= N_x - k - 1; --i) {
				Keep(i, j, false, num_keep);
			}
		}

		if (num_keep[0] < max_keep) {
			m = (max_keep - num_keep[0]) / 4;
			n = 0;
			for (n = 0; n < m; ++n) {
				k = uniformRandom(N);
//				j = invert(k, i);
				i = k / N_y; // NOTE INLINED INVERT
				j = k - i * N_y; // NOTE INLINED INVERT
				Keep(i, j, true, num_keep);
			}
		}

		i = uniformRandom(N_x);
		j = uniformRandom(N_y);

		if (num_keep[0] < max_keep) {
			for (ii = i - N_x / 30; ii < (i + N_x / 30); ++ii) {
				for (jj = j - N_y / 30; jj < (j + N_y / 30); ++jj) {
					Keep(ii, jj, true, num_keep);
					if (num_keep[0] >= max_keep) {
						break;
					}
				}
				if (num_keep[0] >= max_keep) {
					break;
				}
			}
		}

		// 3. select seed cell(s)
		m = 7 * N_x / 9;
		n = 7 * N_y / 9;
		mm = N_x / 9;
		nn = N_y / 9;
		do {
			i = uniformRandom(m);
			j = uniformRandom(n);
			i += mm;
			j += nn;
		} while (!IsPossible(i, j));
		Mark(i, j, num_cells);

		// 4. randomly add cells to already selected cells
		SelectRandomCells(new int[] { 1 }, max_cells, num_keep, max_keep);

		// 5. use current polygon as "seed" for a refinement
		while (hierarchy_cntr <= hierarchy) {
			/*                                                                  */
			/* reset grid data */
			/*                                                                  */
			old_full = full;
			old_keep = keep;
			N_x_old = N_x;
			N_y_old = N_y;
			N_x *= 3;
			N_y *= 3;
			N = N_x * N_y;

			initializeGrid();
			num_candidates = 0;
			num_keep_candidates = 0;
			max_keep = (int) (N * keep_percent);
			max_cells = (int) (N * mark_percent);
			num_keep = new int[1]; // =0
			num_cells = new int[1]; // =0

			/*                                                                  */
			/* copy data from coarse grid to fine grid */
			/*                                                                  */
			for (i = 0; i < N_x_old; ++i) {
				for (j = 0; j < N_y_old; ++j) {
					ii = i * 3;
					jj = j * 3;
					if (IsOldFull(i, j, N_x_old, N_y_old, old_full)) {
						Mark(ii + 1, jj + 1, num_cells);
						if (IsOldFull(i - 1, j, N_x_old, N_y_old, old_full) && IsOldFull(i, j - 1, N_x_old, N_y_old, old_full)) {
							Mark(ii, jj, num_cells);
						}
						if (IsOldFull(i - 1, j, N_x_old, N_y_old, old_full) && IsOldFull(i, j + 1, N_x_old, N_y_old, old_full)) {
							Mark(ii, jj + 2, num_cells);
						}
						if (IsOldFull(i + 1, j, N_x_old, N_y_old, old_full) && IsOldFull(i, j + 1, N_x_old, N_y_old, old_full)) {
							Mark(ii + 2, jj + 2, num_cells);
						}
						if (IsOldFull(i + 1, j, N_x_old, N_y_old, old_full) && IsOldFull(i, j - 1, N_x_old, N_y_old, old_full)) {
							Mark(ii + 2, jj, num_cells);
						}
						if (IsOldFull(i - 1, j, N_x_old, N_y_old, old_full)) {
							Mark(ii, jj + 1, num_cells);
						}
						if (IsOldFull(i + 1, j, N_x_old, N_y_old, old_full)) {
							Mark(ii + 2, jj + 1, num_cells);
						}
						if (IsOldFull(i, j - 1, N_x_old, N_y_old, old_full)) {
							Mark(ii + 1, jj, num_cells);
						}
						if (IsOldFull(i, j + 1, N_x_old, N_y_old, old_full)) {
							Mark(ii + 1, jj + 2, num_cells);
						}
					} else {
						Keep(ii + 1, jj + 1, true, num_keep);
						if (!IsOldFull(i - 1, j, N_x_old, N_y_old, old_full) && !IsOldFull(i, j - 1, N_x_old, N_y_old, old_full)) {
							Keep(ii, jj, true, num_keep);
						}
						if (!IsOldFull(i - 1, j, N_x_old, N_y_old, old_full) && !IsOldFull(i, j + 1, N_x_old, N_y_old, old_full)) {
							Keep(ii, jj + 2, true, num_keep);
						}
						if (!IsOldFull(i + 1, j, N_x_old, N_y_old, old_full) && !IsOldFull(i, j + 1, N_x_old, N_y_old, old_full)) {
							Keep(ii + 2, jj + 2, true, num_keep);
						}
						if (!IsOldFull(i + 1, j, N_x_old, N_y_old, old_full) && !IsOldFull(i, j - 1, N_x_old, N_y_old, old_full)) {
							Keep(ii + 2, jj, true, num_keep);
						}
						if (!IsOldFull(i - 1, j, N_x_old, N_y_old, old_full)) {
							Keep(ii, jj + 1, true, num_keep);
						}
						if (!IsOldFull(i + 1, j, N_x_old, N_y_old, old_full)) {
							Keep(ii + 2, jj + 1, true, num_keep);
						}
						if (!IsOldFull(i, j - 1, N_x_old, N_y_old, old_full)) {
							Keep(ii + 1, jj, true, num_keep);
						}
						if (!IsOldFull(i, j + 1, N_x_old, N_y_old, old_full)) {
							Keep(ii + 1, jj + 2, true, num_keep);
						}
					}
				}
			}

//     old_full = FreeBMatrix(old_full, N_x_old);
//     old_keep = FreeBMatrix(old_keep, N_x_old);

			if (hierarchy_cntr <= hierarchy) {
				// printf("new round: num_keep = %d, max_keep = %d, num_cells = %d, max_cells =
				// %d\n", num_keep, max_keep, num_cells, max_cells);
				SelectRandomCells(num_cells, max_cells, num_keep, max_keep);
			} else {
				// printf("thinned: num_keep = %d, max_keep = %d, num_cells = %d, max_cells =
				// %d\n", num_keep, max_keep, num_cells, max_cells);
			}

			hierarchy_cntr += 1;
		}

		/**************************************************************************/
		/*                                                                        */
		/* output polygon */
		/*                                                                        */
		/**************************************************************************/
		/*                                                                        */
		/* determine boundaries between marked and unmarked cells */
		/*                                                                        */
		edges = MakeEMatrix(N_x + 1, N_y + 1);

		if (diagonal) {
			top = makeBMatrix(N_x, N_y, false);
			bot = makeBMatrix(N_x, N_y, false);
			lft = makeBMatrix(N_x, N_y, false);
			rgt = makeBMatrix(N_x, N_y, false);
			for (i = 0; i < N_x; ++i) {
				for (j = 0; j < N_y; ++j) {
					if (IsFull(i, j)) {
						setBot(i, j, true);
						setTop(i, j, true);
						setLft(i, j, true);
						setRgt(i, j, true);
					}
				}
			}
			for (i = 0; i <= N_x; ++i) {
				for (j = 0; j <= N_y; ++j) {
					if (isCompletelyFull(i, j)) {
						if (IsFull(i - 1, j) && IsFull(i - 1, j + 1) && IsFull(i, j + 1) && isCompletelyEmpty(i, j - 1)
								&& isCompletelyEmpty(i + 1, j)) {
							setBot(i, j, false);
							setRgt(i, j, false);
						} else if (IsFull(i + 1, j) && IsFull(i + 1, j + 1) && IsFull(i, j + 1) && isCompletelyEmpty(i, j - 1)
								&& isCompletelyEmpty(i - 1, j)) {
							setBot(i, j, false);
							setLft(i, j, false);
						} else if (IsFull(i + 1, j) && IsFull(i + 1, j - 1) && IsFull(i, j - 1) && isCompletelyEmpty(i, j + 1)
								&& isCompletelyEmpty(i - 1, j)) {
							setTop(i, j, false);
							setLft(i, j, false);
						} else if (IsFull(i - 1, j) && IsFull(i - 1, j - 1) && IsFull(i, j - 1) && isCompletelyEmpty(i, j + 1)
								&& isCompletelyEmpty(i + 1, j)) {
							setTop(i, j, false);
							setRgt(i, j, false);
						}
					}
				}
			}
			for (i = 0; i <= N_x; ++i) {
				for (j = 0; j <= N_y; ++j) {
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
						if (IsRgtFull(i, j)) {
							storeEdge(i, j + 1, i + 1, j);
						} else if (IsLftFull(i, j)) {
							storeEdge(i, j, i + 1, j + 1);
						}
					} else if (IsBotFull(i, j)) {
						if (IsLftFull(i, j)) {
							storeEdge(i, j + 1, i + 1, j);
						} else if (IsRgtFull(i, j)) {
							storeEdge(i, j, i + 1, j + 1);
						}
					} else {
						System.out.format("cell[%s,%s]: neither top nor bottom is full!\n", i, j);
					}
				}
			}
		} else {
			for (i = 0; i <= N_x; ++i) {
				for (j = 0; j <= N_y; ++j) {
					if (IsFull(i, j)) {
						if (!IsFull(i, j - 1)) {
							storeEdge(i, j, i + 1, j);
						}
						if (!IsFull(i - 1, j)) {
							storeEdge(i, j, i, j + 1);
						}
					} else {
						if (IsFull(i, j - 1)) {
							storeEdge(i, j, i + 1, j);
						}
						if (IsFull(i - 1, j)) {
							storeEdge(i, j, i, j + 1);
						}
					}
				}
			}
		}

		for (i1 = 0; i1 <= N_x; ++i1) {
			for (j1 = 0; j1 <= N_y; ++j1) {
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

		List<List<double[]>> rings = new ArrayList<>();
		loop_cntr = 0;
		for (i1 = 0; i1 <= N_x; ++i1) {
			for (j1 = 0; j1 <= N_y; ++j1) {
				if (edges[i1][j1].i1 != NIL) {
					assert (edges[i1][j1].j1 != NIL);
					assert (edges[i1][j1].i2 != NIL);
					assert (edges[i1][j1].j2 != NIL);
					assert (((edges[edges[i1][j1].i1][edges[i1][j1].j1].i1 == i1) && (edges[edges[i1][j1].i1][edges[i1][j1].j1].j1 == j1))
							|| ((edges[edges[i1][j1].i1][edges[i1][j1].j1].i2 == i1)
									&& (edges[edges[i1][j1].i1][edges[i1][j1].j1].j2 == j1)));
				}
				if (!GetVis(i1, j1)) {
					rings.add(OutputLoop(output, i1, j1, total_number, aligned, perturb, loop_cntr, diagonal, smooth));
					++loop_cntr;
				}
			}
		}

		System.out.format("\n%s edges generated within %d loops\n", total_number[0] - loop_cntr, loop_cntr);
		return rings;
	}

	private int N_x;
	private int N_y;
	private String file_name;
	private boolean perturb;
	private boolean aligned;
	private boolean diagonal;
	private int hierarchy = 0;
	private int smooth = 0;
	private long seed = 0;
	private double mark_percent = 0.5;
	private boolean holes = false;

	private final Random rand;

	Coord[] pnts = new Coord[0];
	int[] candidates = new int[0];
	int[] keep_candidates = new int[0];
	VertexNode[] vertices = new VertexNode[0];

	/**
	 * Generates a random polygon based on a grid with Nx times Ny quadratic cells.
	 * The number of vertices of the polygon generated is random, but it does depend
	 * on Nx, Ny and the percentage mark_percent: The larger Nx and Ny, the more
	 * vertices the polygon tends to have for a given mark_percent.
	 * <p>
	 * Visually pleasing "random" polygons can be achieved by selecting fairly small
	 * values for mark_percent, e.g., mark_percent:=0.1 or even mark_percent:=0.01.
	 * (However, a small value of mark_percent will also reduce the number of
	 * vertices of the polygon.)
	 * 
	 * @param Nx           The number of cells in the X direction of the grid.
	 * @param Ny           The number of cells in the Y direction of the grid.
	 * @param mark_percent The percentage of vertices marked on the grid. The larger
	 *                     the percentage, the more vertices the polygon tends to
	 *                     have.
	 * @param seed         The seed value used to generate the polygon. Two
	 *                     different calls to this method with the same input values
	 *                     will generate the same polygon if they have the same
	 *                     seed.
	 * @param holes        If true, generates a multiply-connected polygonal area.
	 * @param aligned      If true, all vertices will lie on grid points.
	 * @param perturb      If true, the vertices are moved away from the grid points
	 *                     and (most) polygon edges will not be parallel to the
	 *                     coordinate axes.
	 * @param smoothRounds The number of rounds of corner cutting to apply to the
	 *                     polygon generated. A small positive integer value is
	 *                     recommended.
	 * @param hierarchy    The number of rounds of recursive refinement to apply to
	 *                     the polygon generated. A small positive integer value is
	 *                     recommended. This is akin to increasing the depth of
	 *                     fractal curve.
	 * @param diagonal     If true, cuts off some corners by line segments with
	 *                     inclination +/-1, thus generating an octagonal polygon.
	 * @return A random polygon generated based on the input parameters.
	 */
	public SRPolygonGenerator(int Nx, int Ny, double mark_percent, long seed, boolean holes, boolean aligned, boolean perturb, int smoothRounds,
			int hierarchy, boolean diagonal) {

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
		if (seed < 0) {
			seed = 0;
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
//			mark_percent = 0.5;
		}

		this.N_x = Nx;
		this.N_y = Ny;
		this.mark_percent = mark_percent;
		this.seed = seed;
		this.holes = holes;
		this.aligned = aligned;
		this.perturb = perturb;
		this.smooth = smoothRounds;
		this.hierarchy = hierarchy;
		this.diagonal = diagonal;

		rand = new Random(seed);
		if (diagonal && !perturb) {
			this.aligned = true;
		}
		state = (int) seed;
	}

	public List<List<double[]>> getSRPolygon() {
		return Compute(mark_percent, perturb, aligned, hierarchy, diagonal, smooth, file_name);
	}

}

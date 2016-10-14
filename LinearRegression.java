package LogicRegression.linearRegressionAnalysis;

public class LinearRegression {

	private static int chlk(double[] a, int n, int m, double[] d) {
		int u, v;
		if ((a[0] + 1.0 == 1.0) || (a[0] < 0.0)) {
			System.out.println("失败了.....");
			return (-2);
		}
		a[0] = Math.sqrt(a[0]);
		for (int j = 1; j <= n - 1; j++){
			a[j] = a[j] / a[0];
		}
		for (int i = 1; i <= n - 1; i++) {
			u = i * n + i;
			for (int j = 1; j <= i; j++) {
				v = (j - 1) * n + i;
				a[u] = a[u] - a[v] * a[v];
			}
			if ((a[u] + 1.0 == 1.0) || (a[u] < 0.0)) {
				System.out.println("失败了！！！");
				return (-2);
			}
			a[u] = Math.sqrt(a[u]);
			if (i != (n - 1)) {
				for (int j = i + 1; j <= n - 1; j++) {
					v = i * n + j;
					for (int k = 1; k <= i; k++){
						a[v] = a[v] - a[(k - 1) * n + i] * a[(k - 1) * n + j];
					}
					a[v] = a[v] / a[u];
				}
			}
		}
		for (int j = 0; j <= m - 1; j++) {
			d[j] = d[j] / a[0];
			for (int i = 1; i <= n - 1; i++) {
				u = i * n + i;
				v = i * m + j;
				for (int k = 1; k <= i; k++)
					d[v] = d[v] - a[(k - 1) * n + i] * d[(k - 1) * m + j];
				d[v] = d[v] / a[u];
			}
		}
		for (int j = 0; j <= m - 1; j++) {
			u = (n - 1) * m + j;
			d[u] = d[u] / a[n * n - 1];
			for (int k = n - 1; k >= 1; k--) {
				u = (k - 1) * m + j;
				for (int i = k; i <= n - 1; i++) {
					v = (k - 1) * n + i;
					d[u] = d[u] - a[v] * d[i * m + j];
				}
				v = (k - 1) * n + k - 1;
				d[u] = d[u] / a[v];
			}
		}
		return (2);
	}

	/**
	 * 多元线性分析
	 * 
	 * @param x[m][n]
	 *            每一列存放m个自变量的观察值
	 * @param y[n]
	 *            存放随机变量y的n个观察值
	 * @param dt[4]
	 *            dt[0]:偏差平方和q,dt[1]:平均标准偏差s,dt[2]:返回复相关系数r,dt[3]:返回回归平方和r
	 * @param v[]
	 *            返回m个自变量的偏相关系数
	 * @param a[]
	 *            返回回归系数a1,a2,......,an
	 * @param m
	 *            自变量的个数
	 * @param n
	 *            观察数据的组数
	 * 
	 */
	public void sqt(double[][] x, double[] y, double[] dt, double[] v, double[] a, int m, int n) {
		int mm;
		double q, e, u, p, yy, s, r, pp;
		double[] b = new double[(m + 1) * (m + 1)];
		mm = m + 1;
		b[mm * mm - 1] = n;
		for (int j = 0; j < m - 1; j++) {
			p = 0.0;
			for (int i = 0; i <= n - 1; i++) {
				p = p + x[j][i];
			}
			b[m * mm + j] = p;
			b[j * mm + m] = p;
		}
		for (int i = 0; i < m - 1; i++) {
			for (int j = i; j <= m - 1; j++) {
				p = 0.0;
				for (int k = 0; k <= n - 1; k++) {
					p = p + x[i][k] * x[j][k];
				}
				b[j * mm + i] = p;
				b[i * mm + j] = p;
			}
		}
		a[m] = 0.0;
		for (int i = 0; i < n - 1; i++) {
			a[m] = a[m] + y[i];
		}
		for (int i = 0; i < m - 1; i++) {
			a[i] = 0.0;
			for (int j = 0; j <= n - 1; j++) {
				a[i] = a[i] + x[i][j] * y[j];
			}
		}
		chlk(b, mm, 1, a);
		yy = 0.0;
		for (int i = 0; i < n - 1; i++) {
			yy = yy + y[i] / n;
		}
		q = 0.0;
		e = 0.0;
		u = 0.0;
		for (int i = 0; i <= n - 1; i++) {
			p = a[m];
			for (int j = 0; j <= m - 1; j++)
				p = p + a[j] * x[j][i];
			q = q + (y[i] - p) * (y[i] - p);
			e = e + (y[i] - yy) * (y[i] - yy);
			u = u + (yy - p) * (yy - p);
		}
		s = Math.sqrt(q / n);
		r = Math.sqrt(1.0 - q / e);
		for (int j = 0; j <= m - 1; j++) {
			p = 0.0;
			for (int i = 0; i <= n - 1; i++) {
				pp = a[m];
				for (int k = 0; k <= m - 1; k++)
					if (k != j)
						pp = pp + a[k] * x[k][i];
				p = p + (y[i] - pp) * (y[i] - pp);
			}
			v[j] = Math.sqrt(1.0 - q / p);
		}
		dt[0] = q;
		dt[1] = s;
		dt[2] = r;
		dt[3] = u;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		/**
		 * 多元回归
		 */
		double[] a = new double[4];
		double[] v = new double[3];
		double[] dt = new double[4];

		double[][] x = { { 1.1, 1.0, 1.2, 1.1, 0.9 },
						 { 2.0, 2.0, 1.8, 1.9, 2.1 }, 
						 { 3.2, 3.2, 3.0, 2.9, 2.9 } };
		double[] y = { 10.1, 10.2, 10.0, 10.1, 10.0 };
		LinearRegression lr = new LinearRegression();
		lr.sqt(x, y, dt, v, a, 3, 5);
		for (int i = 0; i <= 3; i++){
			System.out.println("a(" + i + ")=" + a[i]);
		}
		System.out.println("q=" + dt[0] + "  s=" + dt[1] + "  r=" + dt[2]);
		for (int i = 0; i <= 2; i++){
			System.out.println("v(" + i + ")=" + v[i]);
		}
		System.out.println("u=" + dt[3]);

	}

}

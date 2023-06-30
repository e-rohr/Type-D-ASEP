static final ForkJoinPool pool = new ForkJoinPool();
volatile double res;
final double matrix[][];
final int level;
public DeterminantTask(final double matrix[][], final int level) {
this.matrix = matrix;
this.level = level;
}
public void setRawResult(Double d) {
this.res = d;
}
public Double getRawResult() {
return res;
}
public boolean exec() {
if (level == 1) { // Trivial case: 1x1 matrix
res = matrix[0][0];
} else if (level == 2) { // Base case: 2x2 matrix
res = matrix[0][0] * matrix[1][1] - matrix[1][0] *
matrix[0][1];
} 
else { // NxN matrix
    res = 0;
    // have a list of tasks to execute at a later time
    final List<ForkJoinTask<Double>> tasks = new LinkedList<ForkJoinTask<Double>>();
    for (int j1 = 0; j1 < level; j1++) {
        final double[][] m = new double[level - 1][];
    for (int k = 0; k < (level - 1); k++)
    m[k] = new double[level - 1];
    for (int i = 1; i < level; i++) {
    int j2 = 0;
    for (int j = 0; j < level; j++) {
    if (j == j1)
    continue;
    m[i - 1][j2] = matrix[i][j];
    j2++;
    }
    }
    final int idx = j1;
    tasks.add(new ForkJoinTask<Double>() {
        double result;
        public Double getRawResult() {
        return result;
        }
        protected void setRawResult(Double value) {
        result = value;
        }
        protected boolean exec() {
        result = Math.pow(-1.0, 1.0 + idx + 1.0) *
        matrix[0][idx] * pool.invoke(new DeterminantTask(m, level - 1));
            return true;
        }
    });
    }
    // once the task is done completing all of the additional tasks, it will add all of the results.
    for (final ForkJoinTask<Double> task :
        invokeAll(tasks)) {
            res += task.getRawResult();
        }
}
    return true;
}
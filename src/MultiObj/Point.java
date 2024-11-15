package MultiObj;

public class Point {
    double[] objectives;

    public Point(double[] objectives) {
        this.objectives = objectives;
    }

    public double getObjective(int index) {
        return objectives[index];
    }
}

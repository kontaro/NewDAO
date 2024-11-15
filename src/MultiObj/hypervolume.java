package MultiObj;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;


public class hypervolume {
	

	public hypervolume(int a){
		
	}

    public  double calculateHypervolume(List<Point> front, double[] referencePoint) {
        // Sort points by the first objective (ascending order)
        front.sort(Comparator.comparingDouble(p -> p.getObjective(0)));

        double hypervolume = 0.0;
        double previousObjective1 = referencePoint[0];

        for (Point point : front) {
            double objective1 = point.getObjective(0);
            double objective2 = point.getObjective(1);
            hypervolume += (previousObjective1 - objective1) * (referencePoint[1] - objective2);
            previousObjective1 = objective1;
        }

        return hypervolume;
    }

    public  void main(String[] args) {
        List<Point> front = new ArrayList<>();
        front.add(new Point(new double[]{1.0, 2.0}));
        front.add(new Point(new double[]{2.0, 1.5}));
        front.add(new Point(new double[]{3.0, 1.0}));

        double[] referencePoint = {4.0, 3.0};

        double hypervolume = calculateHypervolume(front, referencePoint);
        System.out.println("Hypervolume: " + hypervolume);
    }
}





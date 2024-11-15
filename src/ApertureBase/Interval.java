package ApertureBase;

public class Interval {
	int [] interval = new int [2];
	int alpha;
	
	public Interval(int [] interval, int alpha){
		this.interval = interval;
		this.alpha = alpha;
	}

	public int[] getInterval() {
		return interval;
	}

	public void setInterval(int[] interval) {
		this.interval = interval;
	}

	public int getAlpha() {
		return alpha;
	}

	public void setAlpha(int alpha) {
		this.alpha = alpha;
	}
	
	
	
}

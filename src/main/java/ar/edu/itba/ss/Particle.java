package ar.edu.itba.ss;

public class Particle {

    double[] position;
    double[] speed;
    double radius;
    double mass;

    public Particle(double[] position, double radius, double mass) {
        this.position = position;
        this.radius = radius;
        this.mass = mass;
    }
}

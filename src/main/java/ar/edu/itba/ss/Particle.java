package ar.edu.itba.ss;

import java.util.HashSet;
import java.util.Set;

public class Particle implements Cloneable{

    int id;
    double[] position;
    double[] speed;
    double radius;
    double mass;
    double[] acceleration;
    double[] prevAcceleration;
    Integer cell = null;
    Set<Particle> neighbors;
    double pressure;

    public Particle(int id, double[] position, double radius, double mass) {
        this.id = id;
        this.position = position;
        this.radius = radius;
        this.mass = mass;
        this.speed = new double[]{0, 0};
        this.prevAcceleration = new double[]{0, 0};
        this.neighbors = new HashSet<>();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Particle particle = (Particle) o;

        return id == particle.id;
    }

    @Override
    public int hashCode() {
        return id;
    }

    Particle getClone() throws CloneNotSupportedException {
        return (Particle) super.clone();
    }

    public double getDistanceTo(Particle neighbour) {
        return Math.sqrt(Math.pow(this.position[0] - neighbour.position[0], 2) +
                Math.pow(this.position[1] - neighbour.position[1], 2));
    }

    public double[] getRelativeSpeedTo(Particle neighbour) {
        return new double[]{this.speed[0] - neighbour.speed[0],
                this.speed[1] - neighbour.speed[1]};
    }

    public double getSpeedModule() {
        return Math.sqrt(Math.pow(this.speed[0], 2) + Math.pow(this.speed[1], 2));
    }
}

package ar.edu.itba.ss;

import java.util.Set;

public class Beeman implements Integrator{

    private double dt;

    public Beeman(double dt) {
        this.dt = dt;
    }

    @Override
    public void updatePositions(Set<Particle> particles) {

        for (Particle p: particles){

            if (p.acceleration == null){
                p.acceleration = GranularMaterials.forces(p);
            }

            for (int i = 0; i < p.position.length; i++){
                p.position[i] = p.position[i] + p.speed[i] * dt +
                        (2.0 / 3) * p.acceleration[i] * Math.pow(dt, 2) -
                        (1.0 / 6) * p.prevAcceleration[i] * Math.pow(dt, 2);
            }

        }

    }

    @Override
    public void updateSpeeds(Set<Particle> particles){

        for (Particle p : particles) {

            if (GranularMaterials.checkAndResetPosition(p))
                continue;

            double[] oldSpeed = new double[]{p.speed[0], p.speed[1]};

            for (int i = 0; i < p.speed.length; i++){
                p.speed[i] = p.speed[i] + (3.0 / 2) * p.acceleration[i] * dt -
                        (1.0 / 2) * p.prevAcceleration[i] * dt;
            }

            double[] newForce = GranularMaterials.forces(p);

            for (int i = 0; i < p.speed.length; i++){
                p.speed[i] = oldSpeed[i] + (1.0 / 3) * newForce[i] * dt +
                        (5.0 / 6) * p.acceleration[i] * dt -
                        (1.0 / 6) * p.prevAcceleration[i] * dt;
            }

            p.prevAcceleration = p.acceleration;
            p.acceleration = newForce;
        }
    }
}

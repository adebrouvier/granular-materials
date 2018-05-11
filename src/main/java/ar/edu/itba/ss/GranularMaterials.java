package ar.edu.itba.ss;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

import static ar.edu.itba.ss.CliParser.gamma;

public class GranularMaterials {

    private final static double MIN_RADIUS = 0.01;
    private final static double MAX_RADIUS = 0.015;
    private final static double MAX_DIAMETER = MAX_RADIUS*2;
    private final static double MASS = 0.01;
    private final static double KN = Math.pow(10, 5);
    private final static double GRAVITATIONAL_ACCELERATION = 9.8;
    private static CellIndexMethod cellIndexMethod;

    public static void main(String[] args) throws CloneNotSupportedException {
        CliParser.parseOptions(args);
        cellIndexMethod = new CellIndexMethod(CliParser.width, CliParser.height, MAX_RADIUS + 0.001);
        createParticles();
        initSimulation();
    }

    private static void initSimulation(){
        int iterations = 0;
        printParticles(iterations++);

        double dt = 0.01*Math.sqrt(MASS/KN);
        int dt2 = 0;

        cellIndexMethod.setNeighbors();

        for (double t = 0; t < CliParser.time; t+=dt){

            /* Beeman */
            /*TODO: use beeman for speed dependant forces*/
            updatePositions(dt);

            updateSpeeds(dt);

            if (dt2++ % CliParser.dt2 == 0)
                printParticles(iterations++);
        }
    }

    private static void updateSpeeds(double dt) {
        cellIndexMethod.setNeighbors();

        for (Particle p : cellIndexMethod.particles) {

//            if ((p.position[0] > (CliParser.width/2 - CliParser.opening/2) &&
//                    p.position[0] < (CliParser.width/2 + CliParser.opening/2)) &&
//                    p.position[1] < (-CliParser.height/10) ) { /* Reset particles to the top */
//                p.position = newCoords(p.radius);
//                p.speed = new double[2];
//                p.acceleration = new double[2];
//                p.prevAcceleration = new double[2];
//                continue;
//            }
            double[] newForce = forces(p);

            for (int i = 0; i < p.position.length; i++){

                p.speed[i] = p.speed[i] + (1.0 / 3) * newForce[i] * dt +
                        (5.0 / 6) * p.acceleration[i] * dt -
                        (1.0 / 6) * p.prevAcceleration[i] * dt;
            }

            p.prevAcceleration = p.acceleration;
            p.acceleration = newForce;

            cellIndexMethod.putParticle(p);
        }
    }

    private static double[] newCoords(double radius) {

        double newX = radius + Math.random()*(CliParser.width - radius);
        double newY = CliParser.height - Math.random()*MAX_DIAMETER*3;

        return new double[]{newX, newY};
    }


    private static void updatePositions(double dt) {
        for (Particle p: cellIndexMethod.particles){

            if (p.acceleration == null){
                p.acceleration = forces(p);
            }

            for (int i = 0; i < p.position.length; i++){
                p.position[i] = p.position[i] + p.speed[i] * dt +
                        (2.0 / 3) * p.acceleration[i] * Math.pow(dt, 2) -
                        (1.0 / 6) * p.prevAcceleration[i] * Math.pow(dt, 2);
            }

            cellIndexMethod.putParticle(p);

        }
    }

    private static double[] forces(Particle p) {

        double[] force = new double[2];

        /* Lateral Walls */
        /* Left wall */
        if (p.position[0] < p.radius){
            force = lateralWallCollision(p, 0);
        }

        /* Right wall */
        if (p.position[0] > CliParser.width - p.radius){
            double[] newForce = lateralWallCollision(p, CliParser.width);
            force[0] += newForce[0];
            force[1] += newForce[1];
        }

        /*if (p.position[1] < p.radius &&
                (p.position[0] < (CliParser.width/2 - CliParser.opening/2) ||
                        p.position[0] > (CliParser.width/2 + CliParser.opening/2))){*/
        if (p.position[1] < p.radius){

            double superposition = p.radius - p.position[1];

            if (Math.abs(superposition) > 0){

                double dx = 0;
                double dy = p.position[1];

                double mod = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2));
                double ex = (dx/mod);
                double ey = (dy/mod);

                double relativeSpeed = p.speed[0]*ex + p.speed[1]*ey;

                double normalForce = -KN * superposition - gamma * relativeSpeed;

                force[0] += normalForce * ex;
                force[1] += normalForce * ey;

            }
        }

        /* Particle collision */
        for (Particle neighbour : p.neighbors) {

            if (!neighbour.equals(p)){

                double superposition = p.radius + neighbour.radius - p.getDistanceTo(neighbour);

                if (superposition > 0) {

                    double dx = neighbour.position[0] - p.position[0];
                    double dy = neighbour.position[1] - p.position[1];

                    double mod = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2));
                    double ex = (dx / mod);
                    double ey = (dy / mod);

                    double relativeSpeed = (neighbour.speed[0] - p.speed[0]) * ex + (neighbour.speed[1] - p.speed[1]) * ey;

                    double normalForce = -KN * superposition - gamma * relativeSpeed;

                    force[0] += normalForce * ex;
                    force[1] += normalForce * ey;
                }
            }
        }

        force[1] -= p.mass * GRAVITATIONAL_ACCELERATION;

        force[0] = force[0]/p.mass;
        force[1] = force[1]/p.mass;

        return force;
    }

    private static double[] lateralWallCollision(Particle p, double wallCoord) {

        double fx = 0;
        double fy = 0;

        double superposition = p.radius - Math.abs(p.position[0] - wallCoord);
        if (superposition >= 0) {

            double dx = wallCoord - p.position[0];
            double dy = 0;

            double mod = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2));
            double ex = (dx / mod);
            double ey = (dy / mod);

            double relativeSpeed = p.speed[0]*ex + p.speed[1]*ey;
            double normalForce = -KN * superposition - gamma * relativeSpeed;

            fx += normalForce * ex;
            fy += normalForce * ey;
        }

        return new double[]{fx, fy};
    }

    private static void createParticles() {

        double width = CliParser.width;
        double height = CliParser.height;

        double horizontalLimit = Math.floor(width/MAX_DIAMETER);
        double verticalLimit = Math.floor(height/MAX_DIAMETER);

//        int id = 1;
//        for (int i = 0; i < horizontalLimit; i++){
//            for (int j = 0; j < verticalLimit; j++) {
//                double[] position = {MAX_RADIUS + MAX_DIAMETER*i, MAX_RADIUS + MAX_DIAMETER*j};
//                Particle p = new Particle(id++, position, randomRadius(), MASS);
//                cellIndexMethod.putParticle(p);
//            }
//        }

        Particle p = new Particle(1, new double[]{CliParser.width/5, CliParser.height/6}, randomRadius(), MASS);
        cellIndexMethod.putParticle(p);
        p = new Particle(2, new double[]{CliParser.width/5, CliParser.height/4}, randomRadius(), MASS);
        cellIndexMethod.putParticle(p);

    }

    private static double randomRadius() {
        return MIN_RADIUS + Math.random()*(MAX_RADIUS - MIN_RADIUS);
    }

    private static void printParticles(int iteration){
        System.out.println(cellIndexMethod.particles.size());
        System.out.println(iteration);
        for (Particle p: cellIndexMethod.particles)
            System.out.println(p.position[0] + "\t" + p.position[1] + "\t" + p.radius + "\t" + p.speed[0] + "\t" + p.speed[1]);
    }

}

package ar.edu.itba.ss;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import static ar.edu.itba.ss.CliParser.gamma;

public class GranularMaterials {

    private final static double MIN_RADIUS = 0.01;
    private final static double MAX_RADIUS = 0.015;
    private final static double MAX_DIAMETER = MAX_RADIUS*2;
    private final static double MASS = 0.01;
    private final static double KN = Math.pow(10, 4);
    private final static double GRAVITATIONAL_ACCELERATION = 9.8;
    private static double exitParticles;
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

        Integrator integrator = new Beeman(dt);


        PrintWriter flowFile = null;
        try {
            flowFile = new PrintWriter(CliParser.statsFile, "UTF-8");
        } catch (FileNotFoundException | UnsupportedEncodingException e) {
            System.err.println("Could not write to stats file.");
        }

        for (double t = 0; t < CliParser.time; t+=dt){

            integrator.updatePositions(cellIndexMethod.particles);

            updateCells();

            cellIndexMethod.setNeighbors();

            integrator.updateSpeeds(cellIndexMethod.particles);

            updateCells();

            if (dt2++ % CliParser.dt2 == 0) {
                printParticles(iterations++);
                double kineticEnergy = getKineticEnergy();
                flowFile.println(exitParticles/(dt*dt2) + "\t" + kineticEnergy);
                exitParticles = 0;
            }
        }
        flowFile.close();
    }

    private static double getKineticEnergy() {

        double kineticEnergy = 0;

        for (Particle p: cellIndexMethod.particles){
            kineticEnergy += p.getKineticEnergy();
        }

        return kineticEnergy;
    }

    private static void updateCells(){
        for (Particle p: cellIndexMethod.particles) {
            cellIndexMethod.putParticle(p);
        }
    }

    public static boolean checkAndResetPosition(Particle p) {

        if ((p.position[0] > (p.radius + CliParser.width/2 - CliParser.opening/2) &&
                p.position[0] < (CliParser.width/2 + CliParser.opening/2 - p.radius)) &&
                p.position[1] < (-CliParser.height/10) ) { /* Reset particles to the top */
            p.position = newCoords(p.radius);
            p.speed = new double[2];
            p.acceleration = new double[2];
            p.prevAcceleration = new double[2];
            exitParticles++;
            return true;
        }

        return false;
    }

    private static double[] newCoords(double radius) {

        double newX = radius + Math.random()*(CliParser.width - radius);
        double newY = CliParser.height - Math.random()*MAX_DIAMETER*3;

        return new double[]{newX, newY};
    }

    public static double[] forces(Particle p) {

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

        if (p.position[1] < p.radius &&
                (p.position[0] < (p.radius + CliParser.width/2 - CliParser.opening/2) ||
                        p.position[0] > (CliParser.width/2 + CliParser.opening/2 - p.radius))){

            double superposition = p.radius - p.position[1];

            if (Math.abs(superposition) > 0){

                double dx = 0;
                double dy = -Math.abs(p.position[1]);

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

                    double relativeSpeed = (p.speed[0] - neighbour.speed[0]) * ex + (p.speed[1] - neighbour.speed[1]) * ey;

                    double normalForce = -KN * superposition - gamma * relativeSpeed;

                    force[0] += normalForce * ex;
                    force[1] += normalForce * ey;
                }
            }
        }

        p.pressure = Math.sqrt(Math.pow(force[0], 2) + Math.pow(force[1], 2)) /
                (2 * Math.PI * p.radius);
        force[1] -= p.mass * GRAVITATIONAL_ACCELERATION;

        force[0] = force[0]/p.mass;
        force[1] = force[1]/p.mass;

        return force;
    }

    private static double[] lateralWallCollision(Particle p, double wallCoord) {

        double fx = 0;
        double fy = 0;

        double superposition = p.radius - Math.abs(p.position[0] - wallCoord);
        if (superposition > 0) {

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

        int id = 1;
        for (double i = 0; i < verticalLimit; i++){
            for (double j = 0; j < horizontalLimit; j++) {
                double radius = randomRadius();
                double x = radius + MAX_DIAMETER*j + Math.random()*(MAX_DIAMETER - 2*radius);
                double[] position = {x, MAX_RADIUS + MAX_DIAMETER*i};
                Particle p = new Particle(id++, position, radius, MASS);
                cellIndexMethod.putParticle(p);
            }
        }

//        Particle p = new Particle(1, new double[]{MAX_DIAMETER, CliParser.height}, MAX_RADIUS, MASS);
//        p.speed = new double[]{1.0, 0};
//        cellIndexMethod.putParticle(p);
//        p = new Particle(2, new double[]{0.5, CliParser.height}, MAX_RADIUS, MASS);
//        p.speed = new double[]{-1.0, 0};
//        cellIndexMethod.putParticle(p);

    }

    private static double randomRadius() {
        return MIN_RADIUS + Math.random()*(MAX_RADIUS - MIN_RADIUS);
    }

    private static void printParticles(int iteration){
        System.out.println(cellIndexMethod.particles.size());
        System.out.println(iteration);
        for (Particle p: cellIndexMethod.particles)
            System.out.println(p.position[0] + "\t" + p.position[1] + "\t" + p.radius + "\t" + p.pressure);
    }

}

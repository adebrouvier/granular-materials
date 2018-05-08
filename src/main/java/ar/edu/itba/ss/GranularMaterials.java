package ar.edu.itba.ss;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import static ar.edu.itba.ss.CliParser.gamma;

public class GranularMaterials {

    private final static double MIN_RADIUS = 0.01;
    private final static double MAX_RADIUS = 0.015;
    private final static double MAX_DIAMETER = MAX_RADIUS*2;
    private final static double MASS = 0.01;
    private final static double KN = Math.pow(10, 5);
    private final static double KT = 2*KN;
    private final static double GRAVITATIONAL_ACCELERATION = 9.8;
    private static ArrayList<LinkedList<Particle>> cells;
    private static int matrixSize;

    public static void main(String[] args) throws CloneNotSupportedException {
        CliParser.parseOptions(args);
        List<Particle> particles = createParticles();
        initSimulation(particles);
    }

    private static void initSimulation(List<Particle> particles) throws CloneNotSupportedException {
        int iterations = 0;
        printParticles(particles, iterations++);

        double dt = 0.01*Math.sqrt(MASS/KN);
        int dt2 = 0;

        List<Particle> oldParticles = new ArrayList<>();

        for (double t = 0; t < CliParser.time; t+=dt){

            for (Particle p : particles){
                oldParticles.add(p.getClone());
            }

            /* Beeman */
            /*TODO: use beeman for speed dependant forces*/
            updatePositions(dt, particles, oldParticles);

            updateSpeeds(dt, particles);

            if (dt2++ % CliParser.dt2 == 0)
                printParticles(particles, iterations++);

            oldParticles = new ArrayList<>();
        }
    }

    private static void updateSpeeds(double dt, List<Particle> particles) {
        for (Particle p : particles) {

//            if ((p.position[0] > (CliParser.width/2 - CliParser.opening/2) &&
//                    p.position[0] < (CliParser.width/2 + CliParser.opening/2)) &&
//                    p.position[1] < (-CliParser.height/10) ) { /* Reset particles to the top */
//                p.position = newCoords(p.radius);
//                p.speed = new double[2];
//                p.acceleration = new double[2];
//                p.prevAcceleration = new double[2];
//                continue;
//            }

            double[] newForce = forces(p, particles);

            for (int i = 0; i < p.position.length; i++){

                p.speed[i] = p.speed[i] + (1.0 / 3) * newForce[i] * dt +
                        (5.0 / 6) * p.acceleration[i] * dt -
                        (1.0 / 6) * p.prevAcceleration[i] * dt;
            }

            p.prevAcceleration = p.acceleration;

            int newCell = getCellNumber(p);
            if (p.cell != newCell){
                List<Particle> previousList = cells.get(p.cell);
                previousList.remove(p);
                List<Particle> nextList = cells.get(newCell);
                nextList.add(p);
                p.cell = newCell;
            }
        }
    }

    private static double[] newCoords(double radius) {

        double newX = radius + Math.random()*(CliParser.width - radius);
        double newY = CliParser.height - Math.random()*MAX_DIAMETER*3;

        return new double[]{newX, newY};
    }


    private static void updatePositions(double dt, List<Particle> particles, List<Particle> oldParticles) {
        for (Particle p: particles){

            p.acceleration = forces(p, oldParticles);

            for (int i = 0; i < p.position.length; i++){
                p.position[i] = p.position[i] + p.speed[i] * dt +
                        (2.0 / 3) * p.acceleration[i] * Math.pow(dt, 2) -
                        (1.0 / 6) * p.prevAcceleration[i] * Math.pow(dt, 2);
            }

        }
    }

    private static double[] forces(Particle p, List<Particle> particles) {

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
                (p.position[0] < (CliParser.width/2 - CliParser.opening/2) ||
                        p.position[0] > (CliParser.width/2 + CliParser.opening/2))){

            double superposition = p.radius - p.position[1];

            if (superposition >= 0){

                double normalForce = -KN * superposition - gamma * superposition;
                double tangentialForce = 0;//-KT * superposition * p.getSpeedModule();

                double dx = 0;
                double dy = p.position[1];

                double mod = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2));
                double ex = (dx/mod);
                double ey = (dy/mod);

                force[0] += normalForce * ex + tangentialForce * (-ey);
                force[1] += normalForce * ey + tangentialForce * (ex);

            }
        }

        /* Particle collision */
        for (Particle neighbour : particles) {

            if (neighbour.equals(p)){
                break;
            }

            if (p.getDistanceTo(neighbour) < 2*MAX_RADIUS) {

                double superposition = p.radius + neighbour.radius - p.getDistanceTo(neighbour);

                if (superposition < 0)
                    continue;

                double dx = neighbour.position[0] - p.position[0];
                double dy = neighbour.position[1] - p.position[1];

                double mod = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2));
                double ex = (dx/mod);
                double ey = (dy/mod);

                double normalForce = -KN * superposition - gamma * superposition;
                double[] relSpeed = p.getRelativeSpeedTo(neighbour);
                double tangentialForce = 0;//-KT * superposition * (relSpeed[0]*(-ey) + relSpeed[1]*ex);

                force[0] += normalForce * ex + tangentialForce * (-ey);
                force[1] += normalForce * ey + tangentialForce * (ex);
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
            double normalForce = -KN * superposition - gamma * superposition;
            double tangentialForce = 0;//-KT * superposition * p.getSpeedModule();

            double dx = wallCoord - p.position[0];
            double dy = 0;

            double mod = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2));
            double ex = (dx / mod);
            double ey = (dy / mod);

            fx += normalForce * ex + tangentialForce * (-ey);
            fy += normalForce * ey + tangentialForce * (ex);
        }

        return new double[]{fx, fy};
    }

    private static List<Particle> createParticles() {

        double width = CliParser.width;
        double height = CliParser.height;

        double horizontalLimit = Math.floor(width/MAX_DIAMETER);
        double verticalLimit = Math.floor(height/MAX_DIAMETER);

        double minSize = width < (height * 11.0/10.0)  ? width : (height * 11.0/10.0);
        int rowCellsSize = (int) Math.floor(minSize / MAX_RADIUS);
        int cellNumber = rowCellsSize * rowCellsSize;
        matrixSize = rowCellsSize;
        cells = new ArrayList<>(cellNumber);
//        areaLengthX = width / matrixSize;
//        areaLengthY = (height * 11.0/10.0) / matrixSize;

        for (int i = 0; i < cellNumber; i++){
            cells.add(i, new LinkedList<>());
        }

        List<Particle> particles = new ArrayList<>();
        int id = 1;
        for (int i = 0; i < horizontalLimit; i++){
            for (int j = 0; j < verticalLimit; j++) {
                double[] position = {MAX_RADIUS + MAX_DIAMETER*i, MAX_RADIUS + MAX_DIAMETER*j};
                Particle p = new Particle(id++, position, randomRadius(), MASS);
                particles.add(p);
                insertInCell(p);
            }
        }

//        List<Particle> particles = new ArrayList<>();
//        Particle p = new Particle(1, new double[]{CliParser.width/5, CliParser.height}, randomRadius(), MASS);
//        particles.add(p);
//        p = new Particle(2, new double[]{CliParser.width/5, CliParser.height/2}, randomRadius(), MASS);
//        particles.add(p);

        return particles;
    }

    private static double randomRadius() {
        return MIN_RADIUS + Math.random()*(MAX_RADIUS - MIN_RADIUS);
    }

    private static void printParticles(List<Particle> particles, int iteration){
        System.out.println(particles.size());
        System.out.println(iteration);
        for (Particle p: particles)
            System.out.println(p.position[0] + "\t" + p.position[1] + "\t" + p.radius);
    }

    private static void insertInCell(Particle p){
        p.cell = getCellNumber(p);
        List <Particle> cellParticles = cells.get(getCellNumber(p));
        cellParticles.add(p);
    }

    private static int getCellNumber(Particle p){
        double cellX = Math.floor(p.position[0] / (CliParser.width/matrixSize));
        double cellY = Math.floor(p.position[1] + (CliParser.height * 1.0/10.0) / (CliParser.height/matrixSize));
        return (int) (cellY * matrixSize + cellX);
    }
}

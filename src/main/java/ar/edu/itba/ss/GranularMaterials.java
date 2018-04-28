package ar.edu.itba.ss;

import java.util.ArrayList;
import java.util.List;

public class GranularMaterials {

    private final static double MIN_RADIUS = 0.01;
    private final static double MAX_RADIUS = 0.015;
    private final static double MASS = 0.01;
    private final static double KN = Math.pow(10, 5);
    private final static double KT = 2*KN;

    public static void main(String[] args) {
        CliParser.parseOptions(args);
        List<Particle> particles = createParticles();
        initSimulation(particles);
    }

    private static void initSimulation(List<Particle> particles) {
        printParticles(particles);

        double dt = 0.1*Math.sqrt(MASS/KN);

        for (double t = 0; t < CliParser.time; t+=dt){

            for (Particle p: particles){

            }

        }
    }

    private static List<Particle> createParticles() {

        double width = CliParser.width;
        double height = CliParser.height;

        double horizontalLimit = Math.floor(width/(MAX_RADIUS*2));
        double verticalLimit = Math.floor(height/(MAX_RADIUS*2));

        List<Particle> particles = new ArrayList<>();

        for (int i = 0; i < horizontalLimit; i++){
            for (int j = 0; j < verticalLimit; j++) {
                double[] position = {MAX_RADIUS + MAX_RADIUS*2*j, MAX_RADIUS + MAX_RADIUS*2*i};
                Particle p = new Particle(position, randomRadius(), MASS);
                particles.add(p);
            }
        }

        return particles;
    }

    private static double randomRadius() {
        return MIN_RADIUS + Math.random()*(MAX_RADIUS - MIN_RADIUS);
    }

    private static void printParticles(List<Particle> particles){
        System.out.println(particles.size());
        for (Particle p: particles)
            System.out.println(p.position[0] + "\t" + p.position[1] + "\t" + p.radius);
    }
}

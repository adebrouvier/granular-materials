package ar.edu.itba.ss;

import java.util.Set;

public interface Integrator {

    void updatePositions(Set<Particle> particles);

    void updateSpeeds(Set<Particle> particles);

}

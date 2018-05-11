package ar.edu.itba.ss;

import java.util.*;

public class CellIndexMethod{

    private int matrixSizeRows;
    private int matrixSizeColumns;
    private int numberOfCells;
    private double cellLength;
    private HashMap<Integer, List<Particle>> cells = new HashMap<>();
    private final Integer bottomCell = 0;
    private final Integer topCell;
    Set<Particle> particles = new HashSet<>();

    public CellIndexMethod(double width, double height, double radius){
        this.matrixSizeRows = (int) Math.ceil(width / radius);
        this.matrixSizeColumns = (int) Math.ceil(height / radius);
        this.numberOfCells = this.matrixSizeRows * this.matrixSizeColumns;
        this.topCell = this.numberOfCells + 1;
        this.cellLength = radius;
        for (int i = 0; i < numberOfCells + 2; i++){
            cells.put(i, new LinkedList<>());
        }
    }

    private Integer findKeyOfParticle(Particle p){
        int cellX = (int) Math.ceil(p.position[0] / this.cellLength);
        int cellY = (int) Math.ceil(p.position[1] / this.cellLength);
        if (cellX > this.matrixSizeRows || cellX < 0){
            throw new IllegalStateException("Matrix cell wrong.");
        }
        int key = cellY * this.matrixSizeRows + cellX;
        if (key > this.numberOfCells){
            key = this.topCell;
        }
        if (key < 0){
            key = this.bottomCell;
        }
        return key;
    }

    private List<Integer> findNeighborsIndexCells(Particle p){
        List<Integer> neighborCells = new LinkedList<>();
        if (!p.cell.equals(this.bottomCell) && !p.cell.equals(this.topCell)){
            int cellX = p.cell % this.matrixSizeRows;
            int cellY = p.cell / this.matrixSizeRows;
            if (cellX == 0){
                cellX = this.matrixSizeRows;
            }
            if (cellX > this.matrixSizeRows){
                throw new IllegalStateException("Matrix cell wrong.");
            }
            if (cellY + 1 < this.matrixSizeColumns){
                if (cellX + 1 < this.matrixSizeRows){
                    neighborCells.add((cellY+1) * this.matrixSizeRows + cellX + 1);
                }
                if (cellX - 1 > 0){
                    neighborCells.add((cellY+1) * this.matrixSizeRows + cellX - 1);
                }
                neighborCells.add((cellY+1) * this.matrixSizeRows + cellX);
            }else if(cellY + 1 >= this.matrixSizeColumns){
                neighborCells.add(this.topCell);
            }
            if (cellX < this.matrixSizeRows){
                neighborCells.add(cellY * this.matrixSizeRows + cellX + 1);
            }
        }else if(p.cell.equals(this.bottomCell)){
            for (int i = 0; i <= this.matrixSizeRows; i++){
                neighborCells.add(i);
            }
        }

        // No hace falta agregar las de arriba.
        neighborCells.add(p.cell);
        return neighborCells;
    }

    public void putParticle(Particle p){
        particles.add(p);
        if (p.cell != null){
            cells.get(p.cell).remove(p);
        }
        Integer key = this.findKeyOfParticle(p);
        p.cell = key;
        cells.get(key).add(p);
    }

    public void setNeighbors(){
        cells.forEach((k, v) -> {
            for (Particle p: v){
                p.neighbors = new HashSet<>();
            }
            for (Particle p : v){
                List<Integer> index = findNeighborsIndexCells(p);
                List<Particle> neighbors = new LinkedList<>();
                for (Integer i : index){
                    List<Particle> particleList = cells.get(i);
                    neighbors.addAll(particleList);
                    for (Particle particle : particleList){
                        try {
                            p.neighbors.add(particle.getClone());
                        } catch (CloneNotSupportedException e) {
                            e.printStackTrace();
                        }
                    }
                }
                for (Particle n : neighbors){
                    if (!n.equals(p)){
                        try {
                            n.neighbors.add(p.getClone());
                        } catch (CloneNotSupportedException e) {
                            e.printStackTrace();
                        }
                    }
                }
            }
        });
    }

}

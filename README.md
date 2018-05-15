# granular-materials
Granular materials simulation

## Compilation
```
mvn package
```

## Execution

```
time java -jar target/granular-materials-1.0-SNAPSHOT-jar-with-dependencies.jar -w 0.3 -H 0.5 -d 0.165 -t 10 -dt2 1000 -g 20 -sf stats.txt 
```
Parameters:
 * **-d,--opening &lt;arg>**: Size of opening.
 * **-dt2,--deltatime2 &lt;arg>**: Time step for the animation.
 * **-g,--gamma &lt;arg>**: Gammma constant value.
 * **-h, --help**: Shows this screen.
 * **-H,--height &lt;arg>**:Height of the silo.
 * **-sf,--stats &lt;arg>**: Path to the file to output stats.
 * **-t,--time &lt;arg>**: Total time of the simulation.
 * **-w,--width &lt;arg>**: Width of the silo.


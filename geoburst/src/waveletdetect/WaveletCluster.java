package waveletdetect;

import com.mongodb.BasicDBObject;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Created by chao on 9/10/15.
 */
public class WaveletCluster {
    double score;
    List<Integer> entityIds = new ArrayList<Integer>();

    public void setScore(double score) {
        this.score = score;
    }

    public int size() {
        return entityIds.size();
    }

    public void add(int entityId) {
        entityIds.add(entityId);
    }

    public BasicDBObject toBSon() {
        return new BasicDBObject()
                .append("score", score)
                .append("entityIds", entityIds);
    }

    public List<WaveletCluster> split(Map<Integer, Map<Integer, Double>> distanceMap) {
        int randIndex = new Random().nextInt(entityIds.size());
        int entityId = entityIds.get(randIndex);
        int farthestId = entityId;
        for (int iter = 0; iter < 6; iter ++) {
            entityId = farthestId;
            farthestId = findFarthest(distanceMap.get(entityId), entityId);
        }

        WaveletCluster c1 = new WaveletCluster();
        c1.add(entityId);
        WaveletCluster c2 = new WaveletCluster();
        c2.add(farthestId);

        for (int id : entityIds) {
            if (id == entityId || id == farthestId)
                continue;
            assign(c1, c2, entityId, farthestId, id, distanceMap.get(id));
        }
        List<WaveletCluster> ret = new ArrayList<WaveletCluster>();
        ret.add(c1);
        ret.add(c2);
        return ret;
    }

    public int findFarthest(Map<Integer, Double> simMap, int entityId) {
        int farthestId = entityId;
        double minSim = Double.MAX_VALUE;
        for (int other : entityIds) {
            double similarity = simMap.get(other);
            if (similarity < minSim) {
                minSim = similarity;
                farthestId = other;
            }
        }
        return farthestId;
    }

    public void assign(WaveletCluster c1, WaveletCluster c2, int center1, int center2, int id,
                       Map<Integer, Double> distanceMap) {
        double sim1 = distanceMap.get(center1);
        double sim2 = distanceMap.get(center2);
        if (sim1 >= sim2) {
            c1.add(id);
        } else {
            c2.add(id);
        }
    }

}

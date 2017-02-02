package eventweet;

import com.mongodb.BasicDBObject;

import java.util.HashMap;
import java.util.Map;

/**
 * Created by DavidZhou on 7/9/15.
 */
class Descriptors{
    Integer k;
    int o; // # times of bursty and assigned to cluster
    long e; // last bursty frame
    double b_degree;
    double[] spatialDist;

    public Descriptors(Integer keywordID, long fc, double burst_deg, double[] dist)
    {
        k = keywordID;
        o = 1;
        e = fc;
        b_degree = burst_deg;
        spatialDist = dist;
    }

}

public class EventCluster{
    //    List<Entity_Spatial> element;
    double[] centroid; // so called spatial coverage
    long st;
    Map<Integer, Descriptors> Y;
    double score;


    public EventCluster(Integer keywordID, long fc, double b_degree, double[] spatialDist)
    {
//        element = new Vector<Entity_Spatial>();
//        element.add(newKeyword);
        centroid = spatialDist.clone();
        st = fc;

        Y = new HashMap<Integer, Descriptors>();
        Y.put(keywordID, new Descriptors(keywordID, fc, b_degree, spatialDist));
    }

    public void updateCluster(Integer keywordID, long fc, double b_degree, double[] spatialDist)
    {
        if(Y.containsKey(keywordID))
        {
            System.out.println("Same keyword" + keywordID);
            double [] old_dist = Y.get(keywordID).spatialDist;
            for(int i = 0; i < centroid.length; i++)
            {
                centroid[i] = centroid[i]- old_dist[i] + spatialDist[i];
            }
            Y.get(keywordID).b_degree = b_degree;
            Y.get(keywordID).spatialDist = spatialDist;
            Y.get(keywordID).e = fc;
            Y.get(keywordID).o+=1;
        }
        else
        {
            Y.put(keywordID, new Descriptors(keywordID, fc, b_degree, spatialDist));
            int numElem = Y.size() - 1;
            for(int i = 0; i < centroid.length; i++)
            {
                centroid[i] = (centroid[i]*numElem + spatialDist[i])/(numElem + 1);
            }
        }

    }

    public double distance(double[] spatialDist)
    {
        return 1 - cosineSimilarity(this.centroid, spatialDist);
    }
    /* In range [0,1]*/
    private double cosineSimilarity(double[] vectorA, double[] vectorB) {
        double dotProduct = 0.0;
        double normA = 0.0;
        double normB = 0.0;
        for (int i = 0; i < vectorA.length; i++) {
            dotProduct += vectorA[i] * vectorB[i];
            normA += Math.pow(vectorA[i], 2);
            normB += Math.pow(vectorB[i], 2);
        }
        return dotProduct / (Math.sqrt(normA) * Math.sqrt(normB));
    }

    @Override
    public String toString() {
        String s = "Cluster: ";
        s += "size=" + Y.size();
        s += " score=" + score;
        s += " keywords=";
        if(Y.size() > 2) {
            for(Integer entityId: Y.keySet())
                s += " " + entityId.toString();
        }
        return s;
    }


    public BasicDBObject toBSon() {
        return new BasicDBObject()
                .append("score", score)
                .append("size", Y.size())
                .append("entityIds", Y.keySet());
    }

}

package eventweet;

/**
 * Created by DavidZhou on 6/20/15.
 */

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import com.mongodb.BasicDBObject;
import geo.GeoTweet;
import geo.TweetDatabase;
//import sun.org.mozilla.javascript.internal.EcmaError;

/**
 *
 * Created by DavidZhou on 6/13/15.
 */

public class EvenTweet {

    List <TweetDatabase> twTimeGroup = new ArrayList<TweetDatabase>();

    //    private List<Location.> grid = new HashMap<Calendar, List<GeoTweet>>();
    private double[] lngInterval;
    private double[] latInterval;
    private List<List<GeoTweet>> TweetsInGrid;
    double elapseTime;
    //    private List<Entity_Bursty> BurstEntityAndDegree = new ArrayList<Entity_Bursty>();
    private List<EventCluster> ClusterList = new Vector<EventCluster>();

    private Map<Integer, LinkedList<Double>> Hist = new HashMap<Integer, LinkedList<Double>>();
    //    private Map<String, double[]> Hist = new HashMap<String, double[]>();
    private Map<Integer, Double> BurstEntityAndDegree = new HashMap<Integer, Double>();
    private Map<Integer, double[]> BurstEntityAndSpatialDist = new HashMap<Integer, double[]>();
    int numGrid;
    double clusterThreshold, entropyThreshold;

    /* Constructor: read from file to fill up the list of GeoTweet */
    public void selectFromAll(long startTime, long endTime, int windowNum)
    {
        /* Naive way, re assign all tws in current frame */

    }
    public void Init(TweetDatabase refTD, TweetDatabase curTD, int windowNum)
    {
        for (int i = 0; i < windowNum - 1; i++)
        {
            TweetDatabase td = new TweetDatabase();
            twTimeGroup.add(td);
        }
        long windowSize = (refTD.getEndTimestamp() - refTD.getStartTimestamp())/(windowNum-1);
        for (GeoTweet tw : refTD.getTweets())
        {
            int index = (int)((tw.getTimestamp() - refTD.getStartTimestamp())/windowSize);
            if (index >= twTimeGroup.size()) {
                index = twTimeGroup.size() - 1;
            }
            TweetDatabase td = twTimeGroup.get(index);
            td.add(tw);
        }
        twTimeGroup.add(curTD);
        return;
    }

    public EvenTweet(int numGrid, double clusterThreshold, double entropyThreshold)
    {
        this.numGrid = numGrid;
        this.clusterThreshold = clusterThreshold;
        this.entropyThreshold = entropyThreshold; //entropyThreshold = Math.log(numGrid*numGrid/4); // Math.log(2);
    }

    public void detect(TweetDatabase refTD, TweetDatabase curTD)
    {
        long st   = System.currentTimeMillis();

        int windowNum = (int) ((refTD.getEndTimestamp() - refTD.getStartTimestamp()) /
                (curTD.getEndTimestamp() - curTD.getStartTimestamp()));

        Init(refTD, curTD, windowNum);
        /* part for doing burst freq of word */
        Set<Integer> EntityIds = new HashSet<Integer>();
//        tr.Hist.clear();
        int cnt = 0;
        for(GeoTweet tw: curTD.getTweets())
        {
            for(int entityId: tw.getEntityIds())  // for each word of tweets in fc
            {
                cnt++;
                if(!EntityIds.contains(entityId))
                    burstWordNew(windowNum, entityId);
//                else
//                    System.out.println("Dupe" + entityId);
                EntityIds.add(entityId);

            }
        }
//        System.out.println("Count " + cnt + "bursty cnt " + BurstEntityAndDegree.size());


    /* part for doing spatial word  Identification*/
        gridInit();
        int zero_count = 0;
        for(int entityId: BurstEntityAndDegree.keySet())
        {
//            long startTime = System.currentTimeMillis();
            double [] usageDensity = new double[numGrid*numGrid];
            for (int i = 0; i < numGrid; i++)
            {
                for(int j = 0; j < numGrid; j++)
                {
                    usageDensity[i*numGrid + j] = spatialKeywordNew(entityId, i * numGrid + j);
                }

            }

            double entropy  = _compute_entropy(usageDensity);


            if(entropy < entropyThreshold)
                BurstEntityAndSpatialDist.put(entityId, usageDensity);
            else {
                zero_count ++;
//                System.out.println("Pruning......" + entityId + " entropy: " + entropy + " thre: " + entropyThreshold );
            }
        }
//        System.out.println("zero" + zero_count);

//        System.out.println("Spatial Dist");
//        for(String entityId:tr.BurstEntityAndSpatialDist.keySet())
//        {
//            System.out.println(tr.EntityFromId.get(entityId) + " " + tr.BurstEntityAndSpatialDist.get(entityId));
//        }
//
//        System.out.println(tr.BurstEntityAndDegree.size() + "  " + tr.BurstEntityAndSpatialDist.size());

//        System.out.println("........Clustering result........");
        NaiveClustering(curTD.getEndTimestamp());
//        System.out.println("# of clusters after prune " + ClusterList.size());


        long totalTime = System.currentTimeMillis() - st;
        elapseTime = (double)totalTime/1000;
        System.out.println("total running time for eventweet: " + elapseTime + "\n");

        for(int listKey: Hist.keySet())
        {
            LinkedList<Double> list = Hist.get(listKey);
            if(list.size() > 0)
                list.removeFirst();
//            System.out.println(list.size() + " List " + listKey  + " "+ list.toString());
        }

        BurstEntityAndDegree.clear();
        BurstEntityAndSpatialDist.clear();
    }

    /* Initiate the this.twTimeGroup structure here. */

    public void burstWordNew(int windowNum, int entityId)
    {
        /* Algorithm
        * Input: all word in current frame: Wc; Output: subset of Wc: the keywords
        * Total # user in fc: |tweets in twTimeGroup[fc]|
        * Loop though all words w in tweets in fc twTimeGroup.
        * Add to hist_w u(w, fc)
        *
        * */
//        List<Double> curList;
        LinkedList<Double> curList;
        if(!this.Hist.containsKey(entityId))
        {
//            hist = new double[NoGroups];
            LinkedList<Double> hist = new LinkedList<Double>();
            Hist.put(entityId, hist);
//            firstTime = true;
            curList = Hist.get(entityId);
//        if (curList.equals(hist))
//            System.out.println("!!!Equal here");

            for(int groupNo = 0; groupNo < this.twTimeGroup.size(); groupNo++)
            {

                if(this.twTimeGroup.get(groupNo).size() > 0) // has tweet in current time frame
                {
//            System.out.println(i + " " + this.twTimeGroup.get(i).get(0).entity_ID);
                    long numOfUserWithW = 0;
                    long numOfUserTotal = this.twTimeGroup.get(groupNo).size();
                    for (GeoTweet tw : this.twTimeGroup.get(groupNo).getTweets())
                    {
                        for (int id: tw.getEntityIds())
                        {
//                        System.out.println(i + " " + EntityFromId.get(id));
                            if (entityId == id)
                            {
//                        System.out.println(" SAME!!! " + EntityFromId.get(id));
                                numOfUserWithW+=1;
                                break;
                            }
                        }

                    }
                    curList.add((double) numOfUserWithW / numOfUserTotal);
//                System.out.println("Group: " + groupNo + " " +  hist +" " + numOfUserTotal  + " " + numOfUserWithW);
                }

                else
                {
                    curList.add(0.0);
                }
            }
        }
        else
        {
//            System.out.println("Duplicate of old keyword" + entityId + " cur size " + this.Hist.get(entityId).size()) ;

            curList = Hist.get(entityId);

            while(curList.size() < windowNum-1)
            {
                curList.add(0.0); // Add previous missing 0s for re-appearance of keyword
            }
            if(this.twTimeGroup.get(windowNum - 1).size() > 0) // has tweet in current time frame
            {
//            System.out.println(i + " " + this.twTimeGroup.get(i).get(0).entity_ID);
                long numOfUserWithW = 0;
                long numOfUserTotal = this.twTimeGroup.get(windowNum - 1).size();
                for (GeoTweet tw : this.twTimeGroup.get(windowNum - 1).getTweets())
                {
                    for (int id: tw.getEntityIds())
                    {
//                        System.out.println(i + " " + EntityFromId.get(id));
                        if ( entityId == id)
                        {
//                        System.out.println(" SAME!!! " + EntityFromId.get(id));
                            numOfUserWithW+=1;
                            break;
                        }
                    }

                }
                curList.add((double) numOfUserWithW / numOfUserTotal);
//                System.out.println("Group: " + groupNo + " " +  hist +" " + numOfUserTotal  + " " + numOfUserWithW);
            }

            else
            {
                curList.add(0.0);
            }
        }


//        System.out.println(hist);
        double mean, std, b_degree;
        boolean firstTime = true;
        Object[] histArray = curList.toArray();
        if (curList.getLast()!=0.0)
        {
//            System.out.println("last item " + curList.getLast());

            for (int i = 0; i < curList.size()-1; i++) {
                if (!histArray[i].equals(0.0))
                {
//                    System.out.println("array[i] is " +array[i]);
                    firstTime = false;
                }
            }
        }

        if(firstTime)
        {
//            System.out.println("1st time");
            mean = 0;
            std = 0;
        }
        else
        {
            mean = _compute_mean(curList);
            std = _compute_std(curList, mean);
        }

        if(std == 0)
            std = curList.getLast()/ 3;

        b_degree = (curList.getLast() - mean)/std;

        if(b_degree >= 3) {
//        if(!firstTime)
//            System.out.println(entityId + "\t" + "b_degree:" + "\t" + b_degree); //+ " " + mean + " " + std);
//            BurstWords.add(e);
            BurstEntityAndDegree.put(entityId, b_degree);
        }
    }


    public void sortCluster() {
        Comparator<EventCluster> comp  = new Comparator<EventCluster>() {
            public int compare(EventCluster o1, EventCluster o2) {
                if ( o1.score - o2.score < 0)
                    return 1;
                else if (o1.score == o2.score)
                    return 0;
                else
                    return -1;
            }
        };
        Collections.sort(ClusterList,comp);
    }

    public void gridInit()
    {

        /*  Given by Chao   */
//        double maxLng = -73.7;
//        double maxLat = 40.95;
//
//        double minLng = -74.3;
//        double minLat = 40.5;

        double minLng = -74.411;
        double minLat = 40.425;
        double maxLng = -73.598;
        double maxLat = 41.001;



//        System.out.println(minLng + " " + maxLng);
//        System.out.println(minLat + " " + maxLat);

        lngInterval = new double[numGrid+1];
        latInterval = new double[numGrid+1];
        for(int i = 0; i < numGrid+1; i++)
        {
            lngInterval[i] = minLng + (maxLng - minLng)*i/numGrid;
            latInterval[i] = minLat + (maxLat - minLat)*i/numGrid;
//            System.out.print(lngInterval[i] + " ");
        }

        TweetsInGrid = new ArrayList<List<GeoTweet>>(numGrid*numGrid);
        for(int i = 0; i < numGrid*numGrid; i++)
        {
            TweetsInGrid.add(new Vector<GeoTweet>());
        }
//        System.out.println(twTimeGroup.get(twTimeGroup.size() - 1).size());
        for(GeoTweet tw: twTimeGroup.get(twTimeGroup.size() - 1).getTweets())
        {
//            double curLng = tw.geo.getLng(), curLat = tw.geo.getLat();
            int lngIdx = Math.abs(Arrays.binarySearch(lngInterval, tw.getLocation().getLng())) - 2;
            int latIdx = Math.abs(Arrays.binarySearch(latInterval, tw.getLocation().getLat())) - 2;

//            System.out.println(lngIdx + " " + latIdx);
//            if (lngIdx < 0)
//                System.out.println(curLng);
//            if (latIdx < 0)
//                System.out.println(curLat);
            TweetsInGrid.get(lngIdx*numGrid + latIdx).add(tw);
//            System.out.println("ngIdx*numGrid + latIdx" + (lngIdx*numGrid + latIdx));
//            System.out.println(curLng + " " + curLat);
//            System.out.println(lngIdx + " " + latIdx);

        }
    }

    public double spatialKeywordNew(int entityId, int gridIdx)
    {
        // The usage ratio of keyword ki in a cell g ∈ G is the number of users using the keyword in g,
        // compute_entropyd by the total number of users in g.
        double usage_ratio = 0.0;
        int hit_user = 0, tot_user = 0;
        for(GeoTweet tw: TweetsInGrid.get(gridIdx))
        {
            tot_user ++;
            for (int id: tw.getEntityIds())
            {
//                        System.out.println(i + " " + EntityFromId.get(id));
                if (entityId == id)
                {
                    hit_user++;
                }
            }
        }
        if(tot_user!=0)
            usage_ratio = (double)hit_user/tot_user;

//        System.out.println(entityId+ "hit " + hit_user + " tot " + tot_user + " usage_ratio" + usage_ratio);
        return usage_ratio;
    }

    public double ClusterScore(EventCluster c, long fc)
    {
        double[] k_score = new double[c.Y.size()];
        int i = 0;
        double sum = 0;
        for(int entityId: c.Y.keySet())
        {
            Descriptors curY = c.Y.get(entityId);
            double lifeTime = (double) (fc - c.st + 1);
            k_score[i] = curY.o*curY.b_degree*(curY.o/lifeTime) * (1- (fc- curY.e)/lifeTime);
//            k_score[i] = (curY.o/lifeTime) * (1- (fc- curY.e)/lifeTime);
//            System.out.print(k_score[i] + " ");

            sum+=k_score[i];
            i++;
        }

//        System.out.print("score " + sum);
        return sum;
    }

    public void pruneCluster(int topK)
    {
        if (topK < ClusterList.size())
            ClusterList.subList(topK, ClusterList.size()).clear();
    }

    public void NaiveClustering(long fc)
    {

        int count = 0;
        for(int entityId: BurstEntityAndSpatialDist.keySet())
        {
            count ++;

//            if (count % 50 == 0)
//                System.out.println(count);
            double minDist = 1.0;
            int minIndex = 0;
            double[] spatialDist = BurstEntityAndSpatialDist.get(entityId);
            double b_degree = BurstEntityAndDegree.get(entityId);
            if(ClusterList.size() == 0)
                ClusterList.add(new EventCluster(entityId, fc, b_degree, spatialDist));
            else
            {
                for(int i = 0; i < ClusterList.size(); i++)
                {
                    EventCluster c  = ClusterList.get(i);
                    double dist = c.distance(spatialDist);
//                    System.out.println("dist: " + dist);
                    if(dist < minDist)
                    {
                        minDist = dist;
                        minIndex = i;
                    }
                }

                if(minDist < clusterThreshold)
                {
                    ClusterList.get(minIndex).updateCluster(entityId, fc, b_degree, spatialDist);
                }
                else
                {
                    ClusterList.add(new EventCluster(entityId, fc, b_degree, spatialDist));
                }
            }

        }
        for(EventCluster c : ClusterList)
        {
//            System.out.print("size of cluster: " + c.Y.size() + " starttime " + c.st + " => ");
//            for(String entityId: c.Y.keySet())
//                System.out.print(EntityFromId.get(entityId) + "\t");
            c.score = ClusterScore(c, fc);
//            c.score = 0;
//            System.out.println("c score "+c.score);
//            break;
        }
        if (ClusterList.size() == 0)
        {
            System.out.println("No EventCluster");
            return;

        }
//        System.out.println("Num EventCluster " + ClusterList.size());
        System.setProperty("java.util.Arrays.useLegacyMergeSort", "true");
        try {
            sortCluster();
        } catch (Exception e) {
        }

        pruneCluster(100);

//        for(EventCluster c : ClusterList)
//        {
//            System.out.println(c.toString());
//        }


    }

    public void writeEvents(String eventFile) throws Exception {
        BufferedWriter bw = new BufferedWriter(new FileWriter(eventFile, true));
        for(EventCluster c : ClusterList) {
            bw.append(c.toString() + "\n");
        }
        bw.close();
    }

    public void writeStats(String statFile) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(statFile, true));
        bw.write(elapseTime + "\n");
        bw.close();
    }


    // added by ZC.
    public List<BasicDBObject> eventsToBSon() {
        List<BasicDBObject> ret = new ArrayList<BasicDBObject>();
        for (int i = 0; i < ClusterList.size(); i++) {
            EventCluster event = ClusterList.get(i);
            ret.add(event.toBSon());
        }
        return ret;
    }

    public BasicDBObject statsToBSon() {
        return new BasicDBObject().append("time", elapseTime);
    }

    /* Helper functions */
    private static double _compute_mean(List<Double> hist)
    {
        double mean = 0.0;
        for(Iterator<Double> it = hist.iterator(); it.hasNext();)
        {
            mean += it.next();
        }
        return mean/hist.size();

    }
    private static double _compute_std(List<Double> hist, double mean)
    {
        double std = 0.0;
        for(Iterator<Double> it = hist.iterator(); it.hasNext();)
        {
            std+= Math.pow((it.next()-mean), 2);
        }
        return Math.sqrt(std/hist.size());

    }
    private static double _compute_entropy(double[] usageDensity)
    {
        double sum = 0.0, entropy = 0.0;
        for(int i = 0; i < usageDensity.length; i++)
            sum+=usageDensity[i];

//        System.out.println(sum);
        for(int i = 0; i < usageDensity.length; i++)
        {
            usageDensity[i] /= sum;
            if (usageDensity[i]!=0)
                entropy -= usageDensity[i]*Math.log(usageDensity[i]);

            if (usageDensity[i] * Math.log(usageDensity[i]) == Double.NaN)
                System.out.println(usageDensity[i]);
        }

        return entropy;
    }
}


package waveletdetect;

/**
 * Created by DavidZhou on 7/31/15.
 */
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import com.mongodb.BasicDBObject;
import geo.GeoTweet;
import geo.Location;
import geo.TweetDatabase;

//import org.apache.commons.math3.ml.clustering.*;
//import org.apache.commons.math3.ml.distance.DistanceMeasure;
import net.sf.javaml.clustering.*;
import net.sf.javaml.clustering.mcl.DoubleFormat;
import net.sf.javaml.core.*;
import net.sf.javaml.distance.DistanceMeasure;

class LTPair{
    Location l;
    long timestamp;
    public LTPair(Location l, long timestamp)
    {
        this.l = l;
        this.timestamp = timestamp;
    }
}

class XYZ{
    double x,y,z0;
    long z;

    public XYZ(double x, double y, long z)
    {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public XYZ(double x, double y, double z0)
    {
        this.x = x;
        this.y = y;
        this.z0 = z0;
    }

    public void showXYZ()
    {
        System.out.println(x + " " + y + " " + z);
    }
    public void showXYZ0()
    {
        System.out.println(x + " " + y + " " + z0);
    }
}

class Cell{
    int x,y,z;
    public Cell(int x, int y, int z)
    {
        this.x = x;
        this.y = y;
        this.z = z;
    }


}

class Subcomponent{
    List<Cell> cellList;
    XYZ Mean = new XYZ(0.0,0.0,0.0);
    XYZ Std = new XYZ(0.0,0.0,0.0);
    public Subcomponent()
    {
        cellList = new ArrayList<Cell>();
    }
    public void setMeanStd()
    {
        for(Cell c: cellList)
        {
            Mean.x += c.x;
            Mean.y += c.y;
            Mean.z0 += c.z;
        }

        Mean.x /= cellList.size();
        Mean.y /= cellList.size();
        Mean.z0 /= cellList.size();

        for(Cell c: cellList)
        {
            Std.x += (c.x - Mean.x)*(c.x - Mean.x);
            Std.y += (c.y - Mean.y)*(c.y - Mean.y);
            Std.z0 += (c.z - Mean.z0)*(c.z - Mean.z0);
        }
        Std.x /= cellList.size(); Std.x = Math.sqrt(Std.x);
        Std.y /= cellList.size(); Std.y = Math.sqrt(Std.y);
        Std.z0 /= cellList.size(); Std.z0 = Math.sqrt(Std.z0);
    }

    public boolean containsZeroInStd()
    {
        if(Std.x == 0 || Std.y == 0 || Std.z0 == 0)
            return true;
        else
            return false;
    }


}

public class WaveletDetect {
    double intervalX, intervalY;
    long intervalZ;

    Map<Integer, List<LTPair>> Q = new HashMap<Integer, List<LTPair>>();
    Map<Integer, Set<Long>> Q2Ps = new HashMap<Integer, Set<Long>>(); // Mapping from word id->tweets id
    Map<Integer, List<XYZ>> Qxyz = new HashMap<Integer, List<XYZ>>();
    Map<Integer, List<Subcomponent>> Qsc = new HashMap<Integer, List<Subcomponent>>();


    List<WaveletCluster> events = new ArrayList<WaveletCluster>();
    double elapsedTime = 0;

    public WaveletDetect()
    {
//        this.intervalX =0.1;
//        this.intervalY =0.1;
//        this.intervalZ = 86400;
        this.intervalX =0.05;
        this.intervalY =0.05;
        this.intervalZ =5000;
// not good: only 1 cluster
//        this.intervalX =0.01;
//        this.intervalY =0.01;
//        this.intervalZ = 8640;
    }

    public void detect(TweetDatabase refTD)
    {

        for(GeoTweet tw:refTD.getTweets())
        {

            List<Integer> entityIds = tw.getEntityIds();
            LTPair curLT = new LTPair(tw.getLocation(), tw.getTimestamp());
//            System.out.println("Time: " + tw.getTimestamp());
            for(int eid:entityIds)
            {
                Set<Long> tweetIds = Q2Ps.get(eid);
                if(tweetIds == null)
                {
                    tweetIds = new HashSet<Long>();
                    Q2Ps.put(eid, tweetIds);
                }
                tweetIds.add(tw.getTweetId());

                List<LTPair> LTList = Q.get(eid);
                if (LTList == null)
                {
                    LTList = new ArrayList<LTPair>();
                    Q.put(eid, LTList);
                }
                LTList.add(curLT);
            }
        }

        // normalization from Q to Qxyz
        for(Integer q: Q.keySet())
        {
            if(Q.get(q).size()==1)
            {
                XYZ coords = new XYZ(0.0,0.0,0);
                List<XYZ> l = new ArrayList<XYZ>();
                l.add(coords);
                Qxyz.put(q, l);
            }
            else
            {
                double minLat = Q.get(q).get(0).l.getLat();
                double minLng = Q.get(q).get(0).l.getLng();
                long minTimestamp = Q.get(q).get(0).timestamp;

                for(LTPair lt: Q.get(q))
                {
                    if(minLat > lt.l.getLat())
                        minLat = lt.l.getLat();
                    if(minLng > lt.l.getLng())
                        minLng = lt.l.getLng();
                    if(minTimestamp > lt.timestamp)
                        minTimestamp = lt.timestamp;
                }
                List<XYZ> l = new ArrayList<XYZ>();
                for(LTPair lt: Q.get(q))
                {
                    XYZ coords = new XYZ(lt.l.getLat() - minLat,lt.l.getLng() - minLng, lt.timestamp - minTimestamp);
                    l.add(coords);
                }
                Qxyz.put(q, l);
            }

        }
//        System.out.println("Qxyz size: " + Qxyz);

//        System.out.println("Q size: " + Q.size());
        // Quantize the data into cells (compute V)
        for(Integer q: Qxyz.keySet())
        {

            if (Qxyz.get(q).size() > 3)
            {
//                System.out.println("q " + q);
//                System.out.println("q size " + Qxyz.get(q).size());
                double maxX = Qxyz.get(q).get(0).x;
                double maxY = Qxyz.get(q).get(0).y;
                long maxZ = Qxyz.get(q).get(0).z;

                for(XYZ elem: Qxyz.get(q))
                {
//                    elem.showXYZ();
                    if(maxX < elem.x)
                        maxX = elem.x;
                    if(maxY < elem.y)
                        maxY = elem.y;
                    if(maxZ < elem.z)
                        maxZ = elem.z;
                }
//                System.out.format("Max: %d,%d,%d\n", (int)(maxX/intervalX),(int)(maxY/intervalY),(int)(maxZ/intervalZ));
                double[][][] V_C_ijk = new double[(int)(maxX/intervalX) + 1][(int)(maxY/intervalY) + 1][(int)(maxZ/intervalZ) + 1];

                for(XYZ elem: Qxyz.get(q))
                {
//                    elem.showXYZ();
                    V_C_ijk[(int)(elem.x/intervalX)][(int)(elem.y/intervalY)][(int)(elem.z/intervalZ)] += 1;
                }
//                System.out.println("Start wavelet");
                wavelet(V_C_ijk);

                List<Subcomponent> subcomponentList = findConnectedComponentNew(V_C_ijk);
//                System.out.println("subcomponentList size " + subcomponentList.size());
                Qsc.put(q, subcomponentList);

            }


        }

    }

    public void event(){

        final Map<Integer, Map<Integer, Double>> distMapOfMap = new HashMap<Integer, Map<Integer, Double>>();

        List<Integer> oldWordArray = new ArrayList<Integer>(Qsc.keySet());
//        System.out.println("Start event generation... Total tag size: " + oldWordArray.size());
        List<Integer> wordArray = new ArrayList<Integer>();
        for(int i:Qsc.keySet())
        {
            if(Qsc.get(i).size() > 0)
            {
                wordArray.add(i);
            }
        }

        System.out.println("word array size " + wordArray.size());
//        System.out.println("After removing insignificant... Total tag size: " + wordArray.size());
        for(int i = 0; i < wordArray.size() ; i++)
        {
//            System.out.println("i: " + i);
            Map<Integer, Double> distMap = new HashMap<Integer, Double>();
            int wi = wordArray.get(i);

            for(int j = 0; j < wordArray.size(); j++)
            {
                int wj = wordArray.get(j);
                if(i==j) {
                    distMap.put(wj, 1.0);
                    continue;
                }
                double semSim = SemSim(wi, wj) + 0.01;
//                System.out.println(semSim);
                double spadist = SpaDistMean(wi, wj);
//                System.out.println(spadist);
                double dist = semSim/(1+spadist);
                distMap.put(wj, dist);
//                if(dist > 0)
//                    System.out.format("> 0 Dist for %d,%d,%f, %f ,%f\n",wi, wj, semSim, spadist, dist);
            }

            distMapOfMap.put(wi, distMap);
//            if (i % 100 == 0)
//                System.out.println("word array" + i);
        }

        DistanceMeasure dm = new DistanceMeasure() {
            public double measure(Instance instance1, Instance instance2) {
                int p1 = instance1.get(0).intValue();
                int p2 = instance2.get(0).intValue();
//                System.out.format("p1,p2 for %d,%d \n",p1, p2);
                double dist = 1- distMapOfMap.get(p1).get(p2);
//                if(dist < 1)
//                    System.out.format("< 1 Dist for %d,%d,%f\n",p1, p2, dist);
//                System.out.format("Dist for %d,%d,%f, %f \n", wi, wj, semSim, spadist);
                return dist;
            }

            public boolean compare(double v, double v1) {
                return false;
            }

            public double getMinValue() {
                return 0;
            }

            public double getMaxValue() {
                return 0;
            }
        };


//        for (int i : distMapOfMap.keySet()) {
//            System.out.println("Query ID:" + i);
//            Map<Integer, Double> dist = distMapOfMap.get(i);
////            System.out.println(dist);
//            for (Map.Entry<Integer, Double> e : dist.entrySet()) {
//                if (e.getValue() != 0) {
//                    System.out.print(e + " ");
//                }
//            }
//            System.out.println();
//        }

//         clustering
        WaveletCluster initialCluster = new WaveletCluster();
        for (int wordId : wordArray)
            initialCluster.add(wordId);

        List<WaveletCluster> splitResult = new ArrayList<WaveletCluster>();
        List<WaveletCluster> beforeSplit = new ArrayList<WaveletCluster>();
        beforeSplit.add(initialCluster);
        List<WaveletCluster> afterSplit = new ArrayList<WaveletCluster>();

//        int cnt = 0;
        while(beforeSplit.size() > 0) {
//            System.out.println("cnt" + cnt++);
            for (WaveletCluster c : beforeSplit) {
                List<WaveletCluster> children = c.split(distMapOfMap);
                for (WaveletCluster child:children)
                {
                    if (child.size() >= 5 && child.size() <= 30) {
                        splitResult.add(child);
                    }
                    else if (child.size() >= 30)
                    {
                        afterSplit.add(child);
                    }

                }
            }
            beforeSplit = afterSplit;
            afterSplit = new ArrayList<WaveletCluster>();
        }

        events = splitResult;

        for (WaveletCluster e :events) {
            System.out.println(e.entityIds);
        }

//        System.out.println("word array:" + wordArray);
//        DensityBasedSpatialClustering dbscan = new DensityBasedSpatialClustering(0.99, 2, dm);
//        Dataset data = loadDataset(wordArray);
//        System.out.println("dbscan data size " + data.size());
////        System.out.println("data:" + data);
//        Dataset[]  clusters = dbscan.cluster(data);
//        System.out.println("Num of Ds:" + clusters.length);
//        for(Dataset d: clusters) {
////            System.out.println("Num of content in ds:" + d.toString() + d.size());
//            WaveletCluster c = new WaveletCluster();
//            System.out.println("cluster: ");
//            c.setScore(d.size()); // use the number of entities as the score for the cluster.
//            for (int i = 0; i < d.size(); i++) {
//                int entityId = (int) d.instance(i).value(0);
//                System.out.print(entityId + " ");
//                c.add(entityId);
//            }
//            events.add(c);
//            System.out.println("\n");
//        }

    }



    Dataset loadDataset(List<Integer> wordArray)
    {
        Dataset d = new DefaultDataset();
        for(Integer elem: wordArray)
        {
            double[] values = new double[] {(double)(elem)};
            Instance e = new DenseInstance(values);
            d.add(e);
        }
        return d;

    }

    public void detectionMain(TweetDatabase refTD)
    {
        System.out.println("\nStart WaveletDetect......");
        System.out.println("refTD size: " + refTD.size());

        long st = System.currentTimeMillis();
        detect(refTD);
        event();
        long totalTime = System.currentTimeMillis() - st;
        elapsedTime = (double)totalTime/1000;
        System.out.println("total running time for Wavelet detect:" + elapsedTime + "\n");
        System.out.println("Num of events:" + events.size() + "\n");
//        for(WaveletCluster c : events)
//            System.out.println(c.toBSon());
    }

//    public static void main(String[] argv)
//    {
//        System.out.println("here," + (int)(3.1/1.5));
//    }
    private List<Subcomponent> findConnectedComponent(double[][][] arr)
    {
        List<Subcomponent> subcomponentList = new ArrayList<Subcomponent>();
        for(int i = 0; i <  arr.length; i++)
        {
            for(int j = 0; j < arr[0].length; j++)
            {
                for(int k = 0; k < arr[0][0].length; k++)
                {
                    if(arr[i][j][k] != 0)
                    {
                        Subcomponent sc = new Subcomponent();
                        // Current Naive: Treat every non-zero cell as one indep sub-component
                        sc.cellList.add(new Cell(i,j,k));
                        sc.setMeanStd();


                        subcomponentList.add(sc);
//                        System.out.format("ijk: %d,%d,%d,%d \t", i,j,k, arr[i][j][k]);
//                        sc.Std.showXYZ0();
                    }

                }
            }
        }
        return subcomponentList;


    }

    private List<Subcomponent> findConnectedComponentNew(double[][][] arr)
    {
//        int[][][] covered = arr.clone();
        LinkedList<Cell> queue = new LinkedList<Cell>();

        List<Subcomponent> subcomponentList = new ArrayList<Subcomponent>();
        List<Cell> foreground = new ArrayList<Cell>();
        for(int i = 0; i <  arr.length; i++)
        {
            for(int j = 0; j < arr[0].length; j++)
            {
                for(int k = 0; k < arr[0][0].length; k++)
                {
                    if(arr[i][j][k] != 0)
                    {
                        foreground.add(new Cell(i,j,k));
                    }
                }
            }
        }



        for(Cell elem: foreground)
        {
            if(arr[elem.x][elem.y][elem.z] > 0)
            {
                Subcomponent sc = new Subcomponent();
                arr[elem.x][elem.y][elem.z] = 0; // set off the foreground of this pixel
                sc.cellList.add(elem);
                queue.addLast(elem);

                while(!queue.isEmpty())
                {
                    Cell c = queue.poll();
                    List<Cell> neighbors = genNeighbors(c, arr.length, arr[0].length, arr[0][0].length);
                    for (Cell nb: neighbors)
                    {
                        if(arr[nb.x][nb.y][nb.z] > 0)
                        {
//                            System.out.format("ijk: %d,%d,%d,%d \t", nb.x,nb.y ,nb.z, arr[nb.x][nb.y][nb.z]);
                            arr[nb.x][nb.y][nb.z] = 0; // set off the foreground of this pixel
                            sc.cellList.add(nb);
                            queue.addLast(nb);
                        }
                    }
                }
                sc.setMeanStd();

//                System.out.println("M: ");
//                sc.Mean.showXYZ0();
//                System.out.println("STD: ");
//                sc.Std.showXYZ0();

                subcomponentList.add(sc);
//                System.out.println("SSC of this queue: " + sc.cellList.size());
            }
        }

        return subcomponentList;

    }

    private List<Cell> genNeighbors(Cell c, int xMax, int yMax, int zMax)
    {
        List<Cell> neighbors = new ArrayList<Cell>();
        for(int x = Math.max(0, c.x - 1); x < Math.min(xMax, c.x+2); x++)
        {
            for(int y = Math.max(0, c.y - 1); y < Math.min(yMax, c.y+2); y++)
            {
                for(int z = Math.max(0, c.z - 1); z < Math.min(zMax, c.z+2); z++)
                {
                    if(x == c.x && y == c.y && z == c.z)
                    {
                        continue;
                    }
                    neighbors.add(new Cell(x,y,z));
                }
            }
        }

        return neighbors;
    }

    private double SemSim(int qi, int qj)
    {
        Set<Long> si =  Q2Ps.get(qi);
        Set<Long> sj =  Q2Ps.get(qj);
//        System.out.println("Sem sim" + si.size() + " " + sj.size());
        Set<Long> sij = new HashSet<Long>(si);
        sij.retainAll(sj);
//        if(sij.size()!=0)
//        {
//            System.out.format("Sem sim tag: %d, %d; si %d, sj %d, sij %d \n", qi, qj, si.size(), sj.size(), sij.size());
//
//        }
        return (double)(sij.size())/Math.min(si.size(), sj.size());
    }

    private double SpaDist(int qi, int qj)
    {
        List<Subcomponent> Sqi = Qsc.get(qi);
        List<Subcomponent> Sqj = Qsc.get(qj);
//        System.out.println("For id:" + qi + " " + qj);
        System.out.println("Subcomponent list sizes" + Sqi.size() + " " + Sqj.size());
        double spaDist = 0.0;
        for(Subcomponent sk: Sqi)
        {
//            System.out.print("Mean");
//            sk.Mean.showXYZ0();
//            System.out.print("Std");
//            sk.Std.showXYZ0();
            if(sk.containsZeroInStd())
            {
                return Double.POSITIVE_INFINITY;
            }
            double D_min = computeD(sk, Sqj.get(0));
            for(Subcomponent vl: Sqj)
            {
                if(vl.containsZeroInStd())
                {
                    continue;
                }
                double D_cur = computeD(sk, vl);
                if(D_min > D_cur)
                {
                    D_min = D_cur;
                }
            }
            spaDist+=D_min;
        }

        return spaDist;
    }

    private double SpaDistMean(int qi, int qj)
    {
        List<Subcomponent> Sqi = Qsc.get(qi);
        List<Subcomponent> Sqj = Qsc.get(qj);
//        System.out.println("For id:" + qi + " " + qj);
//        System.out.println("Subcomponent list sizes" + Sqi.size() + " " + Sqj.size());
        double spaDist = 0.0;
        for(Subcomponent sk: Sqi)
        {

            double D_min = computeDMean(sk, Sqj.get(0));
            for(Subcomponent vl: Sqj)
            {
                double D_cur = computeDMean(sk, vl);
                if(D_min > D_cur)
                {
                    D_min = D_cur;
                }
            }
            spaDist+=D_min;
        }

        return spaDist;
    }

    private double computeDMean(Subcomponent sk, Subcomponent vl)
    {
        double dist = Math.sqrt(
                 Math.pow((sk.Mean.x - vl.Mean.x),2)
                +Math.pow((sk.Mean.y - vl.Mean.y),2)
                +Math.pow((sk.Mean.z - vl.Mean.z),2));
        return dist;
    }
    private double computeD(Subcomponent sk, Subcomponent vl)
    {
        return Math.max(KL(sk, vl), KL(vl, sk));
    }

    private double KL(Subcomponent s1, Subcomponent s2)
    {
        double retval = KLN(s1.Mean.x, s1.Std.x, s2.Mean.x, s2.Std.x)
                +KLN(s1.Mean.y, s1.Std.y, s2.Mean.y, s2.Std.y)
                +KLN(s1.Mean.z0, s1.Std.z0, s2.Mean.z0, s2.Std.z0);
        return retval;
    }

    private double KLN(double mi, double sdi, double mj, double sdj)
    {
        double retval = 1.0/2*(Math.log(Math.pow(sdj,2)/Math.pow(sdi,2)) + Math.pow(sdi,2)/Math.pow(sdj,2) + Math.pow((mi - mj),2)/Math.pow(sdj, 2) -1);
        if (Double.isNaN(retval))
            return Double.POSITIVE_INFINITY;
        else
            return retval;

    }

    private void printArray(int[][][] arr)
    {
        System.out.println("i,j,k:" + arr.length + " " + arr[0].length + " " + arr[0][0].length );
        for(int i = 0; i <  arr.length; i++)
        {
            for(int j = 0; j < arr[0].length; j++)
            {
                for(int k = 0; k < arr[0][0].length; k++)
                {
                    System.out.format("%d,%d,%d,%d\n", i,j,k, arr[i][j][k]);
                }
            }
        }
    }

    private void wavelet_dim(double[][][] V_C_ijk, int dim)
    {
        
        //================= z =================//
        int xlen = V_C_ijk.length;
        int ylen = V_C_ijk[0].length;
        int zlen = V_C_ijk[0][0].length;
        int[] xyzLen = {xlen,ylen,zlen};
        double [] agg = new double[xyzLen[dim]];
        if (agg.length < 4)
            return;
        for(int i = 0; i < xyzLen[dim]; i++)
        {
            int sum = 0;
            for(int j = 0; j < xyzLen[(dim+1)%3]; j++)
            {
                for(int k = 0; k < xyzLen[(dim+2)%3]; k++)
                {
                    if(dim == 0)
                        sum+=V_C_ijk[i][j][k];
                    else if(dim == 1)
                        sum+=V_C_ijk[k][i][j];
                    else
                        sum+=V_C_ijk[j][k][i];
                }
            }
            agg[i] = sum;
        }
//        System.out.println("Before wavelet:\n ");
//        for(int i = 0; i < agg.length; i++)
//            System.out.println(agg[i] + " ");
        Daub db4 = new Daub();
        db4.daubTrans(agg);
//        System.out.println("After wavelet:\n ");
//        for(int i = 0; i < agg.length; i++)
//            System.out.println(agg[i] + " ");
        double mean_nonzero = 0;
        int non_zero_cnt = 0;
        for(int i = 0; i < agg.length; i++)
        {
            if(agg[i]!=0)
            {
                mean_nonzero+=agg[i];
                non_zero_cnt+=1;
            }
        }
        mean_nonzero/=non_zero_cnt;
//        System.out.println("agg mean: " + mean_nonzero);
        if(non_zero_cnt != 0)
        {

            // average
            for(int i = 0; i < xyzLen[dim]; i++)
            {
//                if(agg[i] < 0.8 * mean_nonzero)
                if(agg[i] < mean_nonzero)
                {
                    for(int j = 0; j < xyzLen[(dim+1)%3]; j++)
                    {
                        for(int k = 0; k < xyzLen[(dim+2)%3]; k++)
                        {
                            if(dim == 0)
                                V_C_ijk[i][j][k] = 0;
                            else if(dim == 1)
                                V_C_ijk[k][i][j] = 0;
                            else
                                V_C_ijk[j][k][i] = 0;
                        }
                    }
                }
            }
        }

    }

    private void wavelet(double[][][] V_C_ijk)
    {
        wavelet_dim(V_C_ijk, 0);
        wavelet_dim(V_C_ijk, 1);
        wavelet_dim(V_C_ijk, 2);
//        waveletY();
//
    }


    // added by ZC.

    public void writeEvents(String eventFile) throws Exception {
        BufferedWriter bw = new BufferedWriter(new FileWriter(eventFile, true));
        for(WaveletCluster c : events) {
            bw.append(c.toString() + "\n");
        }
        bw.close();
    }

    public void writeStats(String statFile) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(statFile, true));
        bw.write(elapsedTime + "\n");
        bw.close();
    }

    public List<BasicDBObject> eventsToBSon() {
        List<BasicDBObject> ret = new ArrayList<BasicDBObject>();
        for (int i = 0; i < events.size(); i++) {
            WaveletCluster event = events.get(i);
            ret.add(event.toBSon());
        }
        return ret;
    }

    public BasicDBObject statsToBSon() {
        return new BasicDBObject().append("time", elapsedTime);
    }

}

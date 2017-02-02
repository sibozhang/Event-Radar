package demo;

import com.mongodb.BasicDBObject;
import com.mongodb.DBObject;
import com.mongodb.util.JSON;
import eventweet.EvenTweet;
import graph.Graph;
import hubseek.Detector;
import waveletdetect.WaveletDetect;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * Created by chao on 11/7/16.
 */
public class IOManager {

    private String expFileName;
    private String queryFileName;
    private String vicinityFileName;

    public IOManager(Map config) throws Exception {
        expFileName = (String) ((Map)config.get("file")).get("experiments");
        queryFileName = (String) ((Map)config.get("file")).get("queries");
        vicinityFileName = (String) ((Map)config.get("file")).get("vicinities");
    }

    // write vicnity information for the entity graph
    public void writeVicinity(Graph graph) throws Exception {
        List<DBObject> docs = new ArrayList<DBObject>();
        for (int nodeId = 0; nodeId < graph.numNode(); nodeId++) {
            BasicDBObject doc = new BasicDBObject().append("ID", new Integer(nodeId).toString());
            List<List<String>> neighbors = new ArrayList<List<String>>();
            for (Map.Entry<Integer, Double> entry : graph.getVicinity(nodeId).entrySet()) {
                List<String> neighbor = new ArrayList<String>();
                neighbor.add(entry.getKey().toString());
                neighbor.add(entry.getValue().toString());
                neighbors.add(neighbor);
            }
            doc.append("Neighbors", neighbors);
            docs.add(doc);
        }
        BufferedWriter bw = new BufferedWriter(new FileWriter(vicinityFileName));
        for(DBObject o : docs) {
            bw.append(JSON.serialize(o) + "\n");
        }
        bw.close();
    }

    public void loadVicinity(Graph graph) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(vicinityFileName));
        Map<Integer, Map<Integer, Double>> vicinity = new HashMap<Integer, Map<Integer, Double>>();
        while(true)  {
            String line = br.readLine();
            if(line == null)    break;
            DBObject obj = (DBObject) JSON.parse(line);
            int nodeId = new Integer((String) obj.get("ID"));
            Map<Integer, Double> neighborMap = new HashMap<Integer, Double>();
            List<List<String>> neighbors = (List<List<String>>) obj.get("Neighbors");
            for (List<String> neighbor : neighbors) {
                int neighborId = new Integer(neighbor.get(0)).intValue();
                double rwr = new Double(neighbor.get(1)).doubleValue();
                neighborMap.put(neighborId, rwr);
            }
            vicinity.put(nodeId, neighborMap);
        }
        graph.setVicinity(vicinity);
    }


    public List<Query> loadBatchQueries(Map config) throws Exception {
        int refWindowSize = (Integer) ((Map) config.get("query")).get("refWindowSize");
        int minSup = (Integer) ((Map) config.get("query")).get("minSup");
        List<Query> queries = new ArrayList<Query>();
        BufferedReader br = new BufferedReader(new FileReader(queryFileName));
        while(true)  {
            String line = br.readLine();
            if(line == null)
                break;
            Scanner sr = new Scanner(line);
            if(!sr.hasNextLine())
                break;
            sr.useDelimiter("\\t");
            int startTS = sr.nextInt(); // node id
            int endTS = sr.nextInt();
            queries.add(new Query(startTS, endTS, refWindowSize, minSup));
        }
        br.close();
        return queries;
    }


    public void writeExp(Detector hubseek, EvenTweet et, WaveletDetect wd) throws Exception {
        // get current time when finished running the experiments
        DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        Date date = new Date();
        BasicDBObject doc = new BasicDBObject().append("time", dateFormat.format(date));
        // hubseek
        BasicDBObject hubseekStats = hubseek.statsToBSon();
        List<BasicDBObject> hubseekEvents = hubseek.eventsToBSon();
        doc.append("hubseek_events", hubseekEvents).append("hubseek_stats", hubseekStats);
        // eventweet
        if (et != null) {
            BasicDBObject eventweetStats = et.statsToBSon();
            List<BasicDBObject> eventweetEvents = et.eventsToBSon();
            doc.append("eventweet_events", eventweetEvents).append("eventweet_stats", eventweetStats);
        }
        // wavelet
        if (wd != null) {
            BasicDBObject waveletStats = wd.statsToBSon();
            List<BasicDBObject> waveletEvents = wd.eventsToBSon();
            doc.append("wavelet_events", waveletEvents).append("wavelet_stats", waveletStats);
        }
        // Write events to file.
        BufferedWriter bw = new BufferedWriter(new FileWriter(expFileName, true));
        bw.append(JSON.serialize(doc) + "\n");
        bw.close();
    }

}


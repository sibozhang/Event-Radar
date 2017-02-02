package demo;

import clustream.Clustream;
import eventweet.EvenTweet;
import geo.GeoTweet;
import geo.TweetDatabase;
import hubseek.BatchDetector;
import hubseek.Detector;
import waveletdetect.WaveletDetect;

import java.util.List;
import java.util.Map;

public class Demo {

    static Database db;
    static Clustream clustream;
    static Map config;
    static IOManager mongo;

    /** ---------------------------------- Initialize ---------------------------------- **/
    static void init(String paraFile) throws Exception {
        config = new Config().load(paraFile);
        mongo = new IOManager(config); // init the connection to mongo db.
        initDatabase();  // load the data
        initClustream(); // initialize the clustream
    }


    // load the initial data, and get the graph ready for computing the similarity between entities.
    static void initDatabase() throws Exception {
        // input and output files
        String tweetFile = (String) (((Map)config.get("file")).get("tweets"));
        String entityFile = (String) (((Map)config.get("file")).get("entities"));
        String entityEdgeFile = (String) (((Map)config.get("file")).get("edges"));
        double epsilon = ((List<Double>) ((Map)config.get("hubseek")).get("epsilon")).get(0);
        int numInitTweets = (Integer) ((Map)config.get("clustream")).get("numInitTweets");
        double pRestart = (Double) ((Map)config.get("clustream")).get("pRestart");
        double errorBound = (Double) ((Map)config.get("clustream")).get("errorBound");
        // load data
        db = new Database();
        db.loadInitialTweets(tweetFile, numInitTweets);
        db.loadEntityGraph(entityFile, entityEdgeFile);
        if ((Boolean) ((Map)config.get("clustream")).get("calcVicinity")) {
            db.getEntityGraph().calcVicinity(epsilon, errorBound, pRestart);
            mongo.writeVicinity(db.getEntityGraph());
        } else {
            mongo.loadVicinity(db.getEntityGraph());
        }
        db.getEntityGraph().printStats();
    }

    // initialize clustream and the entity graph
    static void initClustream() throws Exception {
        // clustream parameters
        int numInitClusters = (Integer) ((Map)config.get("clustream")).get("numInitClusters");
        int numMaxClusters = (Integer) ((Map)config.get("clustream")).get("numMaxClusters");
        int numTweetPeriod = (Integer) ((Map)config.get("clustream")).get("numTweetPeriod");
        int outdatedThreshold = (Integer) ((Map)config.get("clustream")).get("outdatedThreshold");
        // load data
        clustream = new Clustream(numMaxClusters, numTweetPeriod, outdatedThreshold);
        clustream.init(db.getInitialTweets().getTweets(), numInitClusters);
    }


    /** ---------------------------------- run batch ---------------------------------- **/
    static void runBatch() throws Exception {
        TweetDatabase queryDB = new TweetDatabase();
        TweetDatabase refDB = new TweetDatabase();
        List<Query> queries = mongo.loadBatchQueries(config);
        int queryIndex = 0;
        Query query = queries.get(queryIndex);
        GeoTweet tweet;
        while((tweet = db.nextTweet()) != null) {
            addTweet(query, queryDB, refDB, tweet);
            clustream.update(tweet);
            if (tweet.getTimestamp() > query.getEndTS()) {
                System.out.println("Query: " + queryIndex + " startTS: " + query.getStartTS() + " endTS: " + query.getEndTS());
                trigger(query, queryDB, refDB);
                queryIndex ++;
                if (queryIndex < queries.size()) {
                    query = queries.get(queryIndex);
                    queryDB = new TweetDatabase();
                    refDB.deleteFromHead(query.getRefStartTS());
                } else
                    break;
            }
        }
        System.out.println("running batch mode done.");
    }


    static void addTweet(Query query, TweetDatabase queryDB, TweetDatabase refDB, GeoTweet tweet) {
        long ts = tweet.getTimestamp();
        if (ts > query.getStartTS() && ts <= query.getEndTS() && !queryDB.containUserId(tweet))
            queryDB.add(tweet);
        if (ts > query.getRefStartTS() && ts <= query.getRefEndTS())
            refDB.add(tweet);
    }

    static void trigger(Query query, TweetDatabase queryDB, TweetDatabase refDB) throws Exception {
        // skip the empty databases
        if (queryDB.size() == 0)
            return;
        Detector hubseek = runHubSeek(query, queryDB);
        EvenTweet eventweet = runEvenTweet(queryDB, refDB);
        WaveletDetect wavelet = runWavelet(queryDB, refDB);
        writeResults(hubseek, eventweet, wavelet);
        evaluateBandwidth(query, queryDB);
        evaluateEpsilon(query, queryDB);
        evaluateEta(query, queryDB);
    }

    // evaluate hubseek
    static Detector runHubSeek(Query query, TweetDatabase queryDB) {
        Detector detector = new BatchDetector(clustream, db.getEntityGraph());
        if ((Boolean) ((Map)config.get("hubseek")).get("run")) {
            double bandwidth = ((List<Double>) ((Map)config.get("hubseek")).get("bandwidth")).get(0);
            double epsilon = ((List<Double>) ((Map)config.get("hubseek")).get("epsilon")).get(0);
            double eta = ((List<Double>) ((Map)config.get("hubseek")).get("eta")).get(0);
            long refTimeSpan = query.getRefEndTS() - query.getRefStartTS();
            int minSup = query.getMinSup();
            System.out.println("Starting Hubseek");
            detector.detect(queryDB, query.getQueryInterval(), bandwidth, epsilon, minSup, refTimeSpan, eta);
            detector.printStats();
            detector.printEvents();
        }
        return detector;
    }

    // evaluate eventweet.
    static EvenTweet runEvenTweet(TweetDatabase queryDB, TweetDatabase refDB) {
        if (!(Boolean) ((Map)config.get("eventweet")).get("run"))   return null;
        System.out.println("Starting eventweet");
        int numGrid = (Integer) ((Map)config.get("eventweet")).get("numGrid");
        double clusteringThreshold = (Double) ((Map)config.get("eventweet")).get("clusteringThre");
        double entropyThreshold = Math.log(numGrid * numGrid / 4.0);
        EvenTweet evenTweet= new EvenTweet(numGrid, clusteringThreshold, entropyThreshold);
        evenTweet.detect(refDB, queryDB);
        return evenTweet;
    }

    // evaluate wavelet.
    static WaveletDetect runWavelet(TweetDatabase queryDB, TweetDatabase refDB) {
        if (!(Boolean) ((Map)config.get("wavelet")).get("run"))     return null;
        System.out.println("Starting wavelet");
        WaveletDetect waveletDetect = new WaveletDetect();
        waveletDetect.detectionMain(refDB);
        return waveletDetect;
    }

    static void evaluateBandwidth(Query query, TweetDatabase queryDB) throws Exception {
        if ((Boolean) ((Map)config.get("hubseek")).get("evalBandwidth") &&
            (Boolean) ((Map)config.get("hubseek")).get("run")) {
            List<Double> bandwidthList = (List<Double>) ((Map) config.get("hubseek")).get("bandwidth");
            // note that we start from the second value because the default one has already been run.
            for (int i = 1; i < bandwidthList.size(); i++) {
                double bandwidth = bandwidthList.get(i);
                double epsilon = ((List<Double>) ((Map) config.get("hubseek")).get("epsilon")).get(0);
                double eta = ((List<Double>) ((Map) config.get("hubseek")).get("eta")).get(0);
                long refTimeSpan = query.getRefEndTS() - query.getRefStartTS();
                int minSup = query.getMinSup();
                Detector detector = new BatchDetector(clustream, db.getEntityGraph());
                detector.detect(queryDB, query.getQueryInterval(), bandwidth, epsilon, minSup, refTimeSpan, eta);
                writeResults(detector, null, null);
            }
        }
    }

    static void evaluateEpsilon(Query query, TweetDatabase queryDB) throws Exception {
        if ((Boolean) ((Map)config.get("hubseek")).get("evalEpsilon") &&
            (Boolean) ((Map)config.get("hubseek")).get("run")) {
            List<Double> epsilonList = (List<Double>) ((Map) config.get("hubseek")).get("epsilon");
            // note that we start from the second value because the default one has already been run.
            for (int i = 1; i < epsilonList.size(); i++) {
                double bandwidth = ((List<Double>) ((Map) config.get("hubseek")).get("bandwidth")).get(0);
                double epsilon = epsilonList.get(i);
                double eta = ((List<Double>) ((Map) config.get("hubseek")).get("eta")).get(0);
                long refTimeSpan = query.getRefEndTS() - query.getRefStartTS();
                int minSup = query.getMinSup();
                Detector detector = new BatchDetector(clustream, db.getEntityGraph());
                detector.detect(queryDB, query.getQueryInterval(), bandwidth, epsilon, minSup, refTimeSpan, eta);
                writeResults(detector, null, null);
            }
        }
    }


    static void evaluateEta(Query query, TweetDatabase queryDB) throws Exception {
        if ((Boolean) ((Map)config.get("hubseek")).get("evalEta") &&
            (Boolean) ((Map)config.get("hubseek")).get("run")) {
            List<Double> etaList = (List<Double>) ((Map) config.get("hubseek")).get("eta");
            // note that we start from the second value because the default one has already been run.
            for (int i = 1; i < etaList.size(); i++) {
                double bandwidth = ((List<Double>) ((Map) config.get("hubseek")).get("bandwidth")).get(0);
                double epsilon = ((List<Double>) ((Map) config.get("hubseek")).get("epsilon")).get(0);
                double eta = etaList.get(i);
                long refTimeSpan = query.getRefEndTS() - query.getRefStartTS();
                int minSup = query.getMinSup();
                Detector detector = new BatchDetector(clustream, db.getEntityGraph());
                detector.detect(queryDB, query.getQueryInterval(), bandwidth, epsilon, minSup, refTimeSpan, eta);
                writeResults(detector, null, null);
            }
        }
    }


    static void writeResults(Detector hubseek, EvenTweet evenTweet, WaveletDetect wavelet) throws Exception {
        mongo.writeExp(hubseek, evenTweet, wavelet);
    }


    /** ---------------------------------- main ---------------------------------- **/
    public static void main(String [] args) throws Exception {
        String paraFile = args.length > 0 ? args[0] : "../run/sample.yaml";
        init(paraFile);
        runBatch();
    }

}

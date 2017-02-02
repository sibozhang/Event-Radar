package clustream;

import geo.Location;

import java.util.*;

public class Snapshot {

	int order;
	long timeFrameId;
	long timestamp;
	Map<Integer, MicroCluster> clusters; //key: clusterId; value: cluster
	
	public Snapshot(int order, long timeFrameId, long timestamp, Map<Integer, MicroCluster> clusters) {
		this.order = order;
		this.timeFrameId = timeFrameId;
		this.timestamp = timestamp;
		this.clusters = new HashMap<Integer, MicroCluster>();
        for (Map.Entry e : clusters.entrySet()) {
            int clusterId = (Integer) e.getKey();
            MicroCluster cluster = (MicroCluster) e.getValue();
            this.clusters.put(clusterId, new MicroCluster(cluster));
        }
	}

	public Map<Integer, MicroCluster> getClusters() {
        return clusters;
	}

	public long getTimestamp() {
		return timestamp;
	}

	// get the different clusters by subtracting a previous snapshot
	public Set<MicroCluster> getDiffClusters(Snapshot prevSnapshot) {
		Map<Integer, MicroCluster> beforeMap = prevSnapshot.getClusters();
		Map<Integer, MicroCluster> endMap = this.clusters;
		Set<MicroCluster> diffSet = new HashSet<MicroCluster>();
		for (MicroCluster originalCluster : endMap.values()) {
			// copy the original cluster to the base to be subtracted, so that the original cluster remains unchanged.
			MicroCluster base = new MicroCluster(originalCluster);
			if (base.isSingle()) {
				/* 1. if it is a single cluster, we have two cases:
				 * 1) beforeMap does not include this cluster, then this must be a new cluster, we keep all the elements in base.
				 * 2) beforeMap includes this cluster, then do the subtraction.
				 */
				if (beforeMap.containsKey(base.clusterID)) {
					MicroCluster before = beforeMap.get(base.clusterID);
					base.subtract(before);
				}
			} else {
				// 2. composite cluster
				Set<Integer> clusterIDSet = base.idSet;
				clusterIDSet.add(base.clusterID);
				for (Integer cid : clusterIDSet) {
					/* clusterIDSet have four cases:
					 * 1. beforeMap contains cid and it is a single cluster, then we can do subtraction directly
					 * 2. beforeMap contains cid and it is a composite cluster, we can also do subtraction directly
					 * 3. beforeMap does not include cid, and cid is a new cluster, no action is needed.
					 * 4. beforeMap does not include cid, but cid is in the idSet of some composite clusters, then
					 *    no action is needed, because it has already been processed in case 1.
					 */
					if (beforeMap.containsKey(cid)) {
						MicroCluster before = beforeMap.get(cid);
						base.subtract(before);
					}
				}
			}
			if (base.num > 0)
				diffSet.add(base);
		}
		return diffSet;
	}


	// get the language model at a specific location, done by retrieving the nearest cluster and compute the tf-idf dist.
	public Map<Integer, Double> genEntityTfIdfDistribution(Location loc) {
		MicroCluster nearestCluster = getNearestMicroCluster(loc);
//		System.out.println("query location:" + loc.toString());
//		System.out.println("nearest cluster:" + nearestCluster.getCentroidLocation().toString());
		Map<Integer, Double> tfDist = nearestCluster.getEntityTFDistribution();
//		System.out.println("tf distribution:" + tfDist);
		Map<Integer, Double> idfDist = getEntityIdfDistribution(nearestCluster);
//		System.out.println("idf distribution:" + idfDist);
		return multiplyTfIdfDist(tfDist, idfDist);
	}

	private MicroCluster getNearestMicroCluster(Location loc) {
		double smallestDist = Double.MAX_VALUE;
		int retClusterId = -1;
		for (Map.Entry<Integer, MicroCluster> e : this.clusters.entrySet()) {
            int clusterId = e.getKey();
			Location current = e.getValue().getCentroidLocation();
			double dist = loc.calcEuclideanDist(current);
			if (dist < smallestDist) {
				smallestDist = dist;
				retClusterId = clusterId;
			}
		}
//		System.out.println("Smallest distance: " + smallestDist);
		return clusters.get(retClusterId);
	}

	private Map<Integer, Double> getEntityIdfDistribution(MicroCluster nearestCluster) {
		int N = clusters.size(); 	// total number of clusters
		Map<Integer, Double> idfs = new HashMap<Integer, Double>();
		Set<Integer> targetEntityIds = nearestCluster.getEntityIds();
		for (Integer entityId : targetEntityIds) {
			int rawIdf = 0;
			for (MicroCluster c : clusters.values()) {
                if (c.containsEntity(entityId))
                	rawIdf ++;
            }
            idfs.put(entityId, (double) rawIdf);
		}
        for (Integer entityId : idfs.keySet()) {
            double n = idfs.get(entityId);
            double idf = Math.log((N - n + 0.5) / (n + 0.5));
            idfs.put(entityId, idf);
        }
		return idfs;
	}

	// multiply tf and idf; then normalize
	private Map<Integer, Double> multiplyTfIdfDist(Map<Integer, Double> tfDist, Map<Integer, Double> idfDist) {
		Map<Integer, Double> tfIdfDist = new HashMap<Integer, Double>();
		double totalWeight = 0; // used for normalization
        for (Map.Entry<Integer, Double> e : tfDist.entrySet()) {
			int entityId = e.getKey();
			double tf = e.getValue();
			double idf = idfDist.get(entityId);
			tfIdfDist.put(entityId, tf * idf);
			totalWeight += tf * idf;
		}
		for (Map.Entry<Integer, Double> e : tfIdfDist.entrySet()) {
			int entityId = e.getKey();
			double tfIdf = e.getValue() / totalWeight;
			tfIdfDist.put(entityId, tfIdf);
		}
		return tfIdfDist;
	}


	@Override
	public String toString() {
		String itemSep = "=";
		StringBuilder sb = new StringBuilder();
		sb.append(Long.toString(this.timestamp));
		sb.append(itemSep + Integer.toString(this.order));
		for (MicroCluster cluster : this.clusters.values())
			sb.append(itemSep + cluster.toString());
		return sb.toString();
	}

}

//	public static Snapshot string2Snapshot(String ss) {
//		String[] items = ss.split(itemSep);
//		long time = new Long(items[0]);
//		int order = new Integer(items[1]);
//		Map<Integer, TweetCluster> map = new HashMap<Integer, TweetCluster>();
//		for (int i = 2; i < items.length; i++) {
//			TweetCluster geoTweetCluster = TweetCluster.string2TCV(items[i]);
//			map.put(geoTweetCluster.clusterID, geoTweetCluster);
//		}
//		return new Snapshot(order, time, map);
//	}

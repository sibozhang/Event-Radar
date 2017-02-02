package eventweet;

import geo.Location;

import java.util.Calendar;
import java.util.List;
import java.util.Vector;

public class Tweets {
    public String tid;
    public String uid;
    public long timestamp;

    public Location geo;
    public List<String> entity_ID = new Vector<String>();


}

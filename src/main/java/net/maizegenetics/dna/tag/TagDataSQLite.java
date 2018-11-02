package net.maizegenetics.dna.tag;

import cern.colt.list.FloatArrayList;
import com.google.common.collect.*;
import com.google.common.io.CharStreams;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.dna.snp.SimpleAllele;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportBuilder;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.db.SQL;
import org.sqlite.SQLiteConfig;

import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Defines xxxx
 * TAS-480
 * @author Ed Buckler
 */
public class TagDataSQLite implements TagDataWriter, AutoCloseable {
    private Connection connection = null;

    /*These maps contain  objects that are most queried by users.  This is a not the simplest way to do this, which
    would probably be done cleaner just with queries against the databases.  However, there are large performance
    gains in the case of SQLite (or at least with my ability to optimize).

    The logic behind this most of the datasets are relatively small (except tagTagIDMap), and this prevents creation
    of these objects over and over again.
     */
    private BiMap<Tag,Integer> tagTagIDMap;
    private BiMap<String,Integer> tissueTissueIDMap;
    private Map<String,Integer> mappingApproachToIDMap;
    private SortedMap<Position,Integer> cutPosToIDMap;
    public BiMap<Position,Integer> snpPosToIDMap;
    private BiMap<Allele,Integer> alleleToIDMap;

    private TaxaList myTaxaList;

    PreparedStatement tagTaxaDistPS;
    PreparedStatement tagAlleleWhereTagPS;
    PreparedStatement tagidWhereSNPidPS;
    PreparedStatement tagidWhereAlleleidPS;
    PreparedStatement posTagInsertPS;
    PreparedStatement taxaDistWhereCutPositionIDPS;
    PreparedStatement snpPositionsForChromosomePS;
    PreparedStatement alleleTaxaDistForSnpidPS;
    PreparedStatement allAlleleTaxaDistForSnpidPS;
    PreparedStatement snpQualityInsertPS;

    public TagDataSQLite(String filename) {
        try{
            Class.forName("org.sqlite.JDBC");
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
            System.err.println(e.getMessage());
        }
        // create a database connection

        try {
            boolean doesDBExist= Files.exists(Paths.get(filename));

            SQLiteConfig config=new SQLiteConfig();
            //config.setSynchronous(SQLiteConfig.SynchronousMode.OFF);
            /*Optimization ideas
            sqlite3_exec(mDb, "PRAGMA synchronous=OFF", NULL, NULL, &errorMessage);
            sqlite3_exec(mDb, "PRAGMA count_changes=OFF", NULL, NULL, &errorMessage);
            sqlite3_exec(mDb, "PRAGMA journal_mode=MEMORY", NULL, NULL, &errorMessage);
            sqlite3_exec(mDb, "PRAGMA temp_store=MEMORY", NULL, NULL, &errorMessage);
             */
            //config.setTempStore(SQLiteConfig.TempStore.MEMORY);
            connection = DriverManager.getConnection("jdbc:sqlite:"+filename,config.toProperties());
            connection.setAutoCommit(true);  //This has massive performance effects
            Statement statement = connection.createStatement();
            statement.setQueryTimeout(30);  // set timeout to 30 sec.
//            System.out.println(schema);
            if(doesDBExist==false) {
                String schema = CharStreams.toString(new InputStreamReader(TagDataSQLite.class.getResourceAsStream("TagSchema.sql")));
                statement.executeUpdate(schema);
            }
            initPreparedStatements();
            loadTagHash();
            loadTissueHash();
            loadMappingApproachHash();
            loadTaxaList();
        }
        catch(Exception e)
        {
            // if the error message is "out of memory",
            // it probably means no database file is found
            System.err.println(e.getMessage());
            e.printStackTrace();
        }
    }

    @Override
    public void close() throws Exception {
        System.out.println("Closing SQLDB");
        connection.close();
    }

    private void initPreparedStatements() {
        try{
            posTagInsertPS=connection.prepareStatement(
                    "INSERT OR IGNORE into tagCutPosition (tagid, positionid, mapappid, bestmapping, forward, cigar, supportval)" +
                            " values(?,?,?,?,?,?,?)");
            tagTaxaDistPS=connection.prepareStatement("select depthsRLE from tagtaxadistribution where tagid=?");
            tagAlleleWhereTagPS=connection.prepareStatement("select * from tagallele where tagid=?");
            tagidWhereSNPidPS=connection.prepareStatement(
                    "select tagid from allele, tagallele where allele.snpid=? and allele.alleleid=tagallele.alleleid");
            tagidWhereAlleleidPS=connection.prepareStatement(
                    "select tagid from tagallele where alleleid=?");
            taxaDistWhereCutPositionIDPS=connection.prepareStatement(
                    "select tagtaxadistribution.* from tagCutPosition, tagtaxadistribution where tagCutPosition.positionid=? and " +
                            "tagCutPosition.tagid=tagtaxadistribution.tagid and tagCutPosition.bestmapping=1");
            snpPositionsForChromosomePS=connection.prepareStatement(
            		"select position, qualityScore, refAllele from snpposition where chromosome=?");
            alleleTaxaDistForSnpidPS =connection.prepareStatement("select a.*, td.* from allele a, tagallele ta, tagtaxadistribution td\n" +
                    "where a.alleleid=ta.alleleid and ta.tagid=td.tagid and a.snpid=?");
            allAlleleTaxaDistForSnpidPS =connection.prepareStatement("select a.*, td.* from allele a, tagallele ta, tagtaxadistribution td\n" +
                    "where a.alleleid=ta.alleleid and ta.tagid=td.tagid order by a.snpid");
            snpQualityInsertPS=connection.prepareStatement(
                    "INSERT into snpQuality (snpid, taxasubset ,avgDepth, minorDepthProp, minor2DepthProp, gapDepthProp, " +
                            "propCovered, propCovered2, taxaCntWithMinorAlleleGE2, minorAlleleFreqGE2, inbredF_DGE2)" +
                            " values(?,?,?,?,?,?,?,?,?,?,?)");
//            snpQualityInsertPS=connection.prepareStatement(
//                    "INSERT into snpQuality (snpid, taxasubset)" +
//                            " values(?,?)");
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void loadTagHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from tag");
            int size=rs.getInt(1);
            System.out.println("size of all tags in tag table=" + size);
            if(tagTagIDMap==null || size/(tagTagIDMap.size()+1)>3) tagTagIDMap=HashBiMap.create(size);
            rs=connection.createStatement().executeQuery("select * from tag");
            boolean hasName;
            try{rs.findColumn("tagName");hasName=true;}catch (SQLException e) {hasName=false;}
            while(rs.next()) {
                TagBuilder tagBuilder=TagBuilder.instance(rs.getBytes("sequence"),rs.getShort("seqlen"));
                if(hasName) tagBuilder.name(rs.getString("tagName"));
                tagTagIDMap.putIfAbsent(tagBuilder.build(),rs.getInt("tagid"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    
    private void loadTissueHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from tissue");
            int size=rs.getInt(1);
            System.out.println("size of all tissues in tissue table=" + size);
            if(tissueTissueIDMap==null || size/(tissueTissueIDMap.size()+1)>3) tissueTissueIDMap=HashBiMap.create(size);
            rs=connection.createStatement().executeQuery("select * from tissue");
            while(rs.next()) {
                tissueTissueIDMap.putIfAbsent(rs.getString("tissue"),rs.getInt("tissueid"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    private void loadCutPositionHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from cutPosition");
            int size=rs.getInt(1);
            System.out.println("size of all positions in cutPosition table="+size);
            if(cutPosToIDMap==null) {cutPosToIDMap=new TreeMap<>();}
            else if(size==cutPosToIDMap.size()) return;
            rs=connection.createStatement().executeQuery("select * from cutPosition");
            while(rs.next()) {
                Position p=new GeneralPosition
                        .Builder(new Chromosome(rs.getString("chromosome")),rs.getInt("position"))
                        .strand(rs.getByte("strand"))
                        .build();
                cutPosToIDMap.putIfAbsent(p, rs.getInt("positionid"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    
    private void loadSNPPositionHash(boolean force) {
    	if (force) { // reload as quality scores have changed.
    		if (snpPosToIDMap != null) snpPosToIDMap.clear(); 
    	}
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from snpposition");
            int size=rs.getInt(1);
            System.out.println("size of all positions in snpPosition table="+size);
            if(snpPosToIDMap==null) {snpPosToIDMap=HashBiMap.create(size);}
            else if(size==snpPosToIDMap.size()) return;
            rs=connection.createStatement().executeQuery("select * from snpposition");
            while(rs.next()) {
                // LCJ - this needs to add .allele(WHICH_ALLELE.Reference,refAllele)
                // WHich means you need a way to know what is the ref allele
                byte refAllele = (byte)rs.getInt("refAllele");
                Position p=new GeneralPosition
                        .Builder(new Chromosome(rs.getString("chromosome")),rs.getInt("position"))
                        .strand(rs.getByte("strand"))
                        .addAnno("QualityScore", rs.getFloat("qualityScore"))
                        .allele(WHICH_ALLELE.Reference,refAllele)
                        .build();
                snpPosToIDMap.putIfAbsent(p, rs.getInt("snpid"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void loadAlleleHash() {
        try{
            loadSNPPositionHash(false);
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from allele");
            int size=rs.getInt(1);
            System.out.println("size of all alleles in allele table="+size);
            if(alleleToIDMap==null) {alleleToIDMap=HashBiMap.create(size);}
            if(size==alleleToIDMap.size()) return;
            rs=connection.createStatement().executeQuery("select * from allele");
            while(rs.next()) {
                int snpid=rs.getInt("snpid");
                Position p=snpPosToIDMap.inverse().get(snpid);
                Allele a=new SimpleAllele(rs.getByte("allelecall"),p);
                alleleToIDMap.putIfAbsent(a, rs.getInt("alleleid"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void loadMappingApproachHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from mappingApproach");
            int size=rs.getInt(1);
            System.out.println("size of all tags in mappingApproach table="+size);
            if(size==0) {
                connection.createStatement().executeUpdate("insert into mappingApproach (approach, software, protocols) " +
                        "values('unknown','unknown','unknown')");
                size=1;
            }
            mappingApproachToIDMap=new HashMap<>(size);
            rs=connection.createStatement().executeQuery("select * from mappingApproach");
            while(rs.next()) {
                mappingApproachToIDMap.put(rs.getString("approach"), rs.getInt("mapappid"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void loadTaxaList() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from taxa");
            int size=rs.getInt(1);
            System.out.println("size of all taxa in taxa table="+size);
            TaxaListBuilder tlb=new TaxaListBuilder();
            rs=connection.createStatement().executeQuery("select * from taxa");
            while(rs.next()) {
                tlb.add(new Taxon(rs.getString("name")));
            }
            myTaxaList=tlb.build();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public TaxaList getTaxaList() {
        if(myTaxaList==null) loadTaxaList();
        return myTaxaList;
    }

    @Override
    public Map<Tag, String> getTagsNameMap() {
        try{
            Map<Tag, String> tagNameMap=new HashMap<>(tagTagIDMap.size()+1);
            ResultSet rs=connection.createStatement().executeQuery("select * from tag");
            while(rs.next()) {
                tagNameMap.put(TagBuilder.instance(rs.getBytes("sequence"), rs.getShort("seqlen")).build(), rs.getString("tagName"));
            }
            return tagNameMap;
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return null;
    }

    @Override
    public boolean putAllTag(Set<Tag> tags) {
        int batchCount=0, totalCount=0;
        try {
            connection.setAutoCommit(false);
            PreparedStatement tagInsertPS=connection.prepareStatement("insert into tag (sequence, seqlen) values(?,?)");
            for (Tag tag : tags) {
                if(tagTagIDMap.containsKey(tag)) continue;  //it is already in the DB skip
                tagInsertPS.setBytes(1, tag.seq2BitAsBytes());
                tagInsertPS.setShort(2, tag.seqLength());
                tagInsertPS.addBatch();
                batchCount++;
                totalCount++;
                if(batchCount>100000) {
                    System.out.println("tagInsertPS.executeBatch() "+batchCount);
                    tagInsertPS.executeBatch();
                    //connection.commit();
                    batchCount=0;
                }
            }
            tagInsertPS.executeBatch();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
            return false;
        }
        if(totalCount>0) loadTagHash();
        return true;
    }

    @Override
    public boolean putAllNamesTag(Map<Tag, String> tagNameMap) {
        int batchCount=0, totalCount=0;
        try {
            connection.setAutoCommit(false);
            PreparedStatement tagInsertPS=connection.prepareStatement("insert into tag (sequence, seqlen, tagName) values(?,?,?)");
            for (Map.Entry<Tag, String> entry : tagNameMap.entrySet()) {
                Tag tag=entry.getKey();
                if(tagTagIDMap.containsKey(tag)) continue;  //it is already in the DB skip
                tagInsertPS.setBytes(1, tag.seq2BitAsBytes());
                tagInsertPS.setShort(2, tag.seqLength());
                tagInsertPS.setString(3, entry.getValue());
                tagInsertPS.addBatch();
                batchCount++;
                totalCount++;
                if(batchCount>100000) {
                    System.out.println("tagInsertPS.executeBatch() "+batchCount);
                    tagInsertPS.executeBatch();
                    //connection.commit();
                    batchCount=0;
                }
            }
            tagInsertPS.executeBatch();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
            return false;
        }
        if(totalCount>0) loadTagHash();
        return true;
    }

    @Override
    public void putTaxaList(TaxaList taxaList) {
        try {
            connection.createStatement().execute("delete from taxa");
            connection.setAutoCommit(false);
            PreparedStatement taxaInsertPS=connection.prepareStatement("insert into taxa (taxonid, name) values(?,?)");
            for (int i = 0; i < taxaList.size(); i++) {
                taxaInsertPS.setInt(1, i);
                taxaInsertPS.setString(2, taxaList.get(i).getName());
                taxaInsertPS.addBatch();

            }
            taxaInsertPS.executeBatch();
            connection.setAutoCommit(true);
            loadTaxaList();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public boolean putTaxaTissueDistribution(String taxon, String tissue, List<Tag> tags, int[] counts) {
        int batchCount=0;
        try {
            connection.setAutoCommit(false);
            int tissueID=tissueTissueIDMap.get(tissue);
            int taxaID=myTaxaList.indexOf(taxon);
            PreparedStatement tagTaxaTissuePS=connection.prepareStatement("replace into tagTaxaTissueDist (tagid, tissueid, taxonid, readCount) values(?,?,?,?) ");
            for (int i = 0; i < counts.length; i++) {
                int tagID=tagTagIDMap.get(tags.get(i));
                tagTaxaTissuePS.setInt(1, tagID);
                tagTaxaTissuePS.setInt(2, tissueID);
                tagTaxaTissuePS.setInt(3, taxaID);
                tagTaxaTissuePS.setInt(4, counts[i]);
                tagTaxaTissuePS.addBatch();
                batchCount++;
                if(batchCount>100000) {
                    System.out.println("putTaxaDistribution next"+batchCount);
                    tagTaxaTissuePS.executeBatch();
                    //connection.commit();
                    batchCount=0;
                }
            }
            tagTaxaTissuePS.executeBatch();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    @Override
    public void putTaxaDistribution(Map<Tag, TaxaDistribution> tagTaxaDistributionMap) {
        int batchCount=0;
        try {
            int numTaxa=myTaxaList.numberOfTaxa();
            connection.setAutoCommit(false);
            PreparedStatement tagInsertPS=connection.prepareStatement("insert into tagtaxadistribution (tagid, depthsRLE, totalDepth) values(?,?,?)");
            for (Map.Entry<Tag, TaxaDistribution> entry : tagTaxaDistributionMap.entrySet()) {
                int tagID=tagTagIDMap.get(entry.getKey());
                tagInsertPS.setInt(1,tagID);
                if(entry.getValue().maxTaxa()!=numTaxa) throw new IllegalStateException("Number of taxa does not agree with taxa distribution");
                tagInsertPS.setBytes(2, entry.getValue().encodeTaxaDepth());
                tagInsertPS.setInt(3, entry.getValue().totalDepth());
                tagInsertPS.addBatch();
                batchCount++;
                if(batchCount>100000) {
                    System.out.println("putTaxaDistribution next"+batchCount);
                    tagInsertPS.executeBatch();
                    //connection.commit();
                    batchCount=0;
                }
            }
            tagInsertPS.executeBatch();
            connection.setAutoCommit(true);  
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void putTagAlignments(Multimap<Tag, Position> tagAnnotatedPositionMap) {
        int batchCount=0;
        try {
            putAllTag(tagAnnotatedPositionMap.keySet());
            putCutPositionsIfAbsent(tagAnnotatedPositionMap.values());
            connection.setAutoCommit(false);
            for (Map.Entry<Tag, Position> entry : tagAnnotatedPositionMap.entries()) {
                Position p=entry.getValue();
                GeneralAnnotation annotation = p.getAnnotation();
                int ind=1;
                posTagInsertPS.setInt(ind++, tagTagIDMap.get(entry.getKey()));
                posTagInsertPS.setInt(ind++, cutPosToIDMap.get(p));
                posTagInsertPS.setInt(ind++, getMappingApproachID(p));
                posTagInsertPS.setBoolean(ind++, true);  //todo this needs to be derived from the position or set later.
                boolean forward=true;
                try{
                    if(annotation.getTextAnnotation("forward")[0].toLowerCase().equals("false")) forward=false;
                } catch (Exception e) {
                    System.err.println(p.toString());
                    System.err.println("Error with forward annotation");
                    //no valid cigarValue
                }
                posTagInsertPS.setBoolean(ind++,forward);
                String cigarValue="";
                try{
                    cigarValue=annotation.getTextAnnotation("cigar")[0];
                } catch (Exception e) {
                    System.err.println(p.toString());
                    System.err.println("Error with cigar");
                    //no valid cigarValue
                }
                posTagInsertPS.setString(ind++, cigarValue);
                short supportVal=0;
                try{
                    String[] svS=annotation.getTextAnnotation("supportvalue");
                    if(svS.length>0) {
                        supportVal=Short.parseShort(svS[0]);
                    }
                } catch (Exception e) {
                    System.err.println("Error with supportVal");
                    //no valid supportVal
                }
                posTagInsertPS.setByte(ind++, (byte) supportVal);
                //System.out.println(posTagInsertPS.toString());
                posTagInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                    System.out.println("putTagAlignments next"+batchCount);
                    posTagInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            posTagInsertPS.executeBatch();
            connection.setAutoCommit(true);
            // print some metrics for debugging
            ResultSet rs = connection.createStatement().executeQuery("select count (DISTINCT positionid) as numCutSites from tagCutPosition");
            if (rs.next()) {
            	System.out.println("Total number of cut sites: " + rs.getInt("numCutSites"));
            }
            PreparedStatement cutSiteNumFromTCPPS = connection.prepareStatement(
                    "select count(*) as numSites from (select count(*) as tgcnt,positionid from tagCutPosition " +
                    "GROUP BY positionid) where tgcnt=?");
            cutSiteNumFromTCPPS.setInt(1, 1);// having 1 tag
            rs = cutSiteNumFromTCPPS.executeQuery();

            if (rs.next()) {
            	System.out.println("Number of cut sites with 1 tag: " + rs.getInt("numSites"));
            }
            cutSiteNumFromTCPPS.setInt(1, 2);// having 2 tag
            rs = cutSiteNumFromTCPPS.executeQuery();
            if (rs.next()) {
            	System.out.println("Number of cut sites with 2 tags: " + rs.getInt("numSites"));
            }
            cutSiteNumFromTCPPS.setInt(1, 3);// having 3 tags
            rs = cutSiteNumFromTCPPS.executeQuery();
            if (rs.next()) {
            	System.out.println("Number of cut sites with 3 tags: " + rs.getInt("numSites"));
            }
            
            PreparedStatement cutSiteGreaterThanPS = connection.prepareStatement(
                    "select count(*) as numSites from (select count(*) as tgcnt,positionid from tagCutPosition " +
                    "GROUP BY positionid) where tgcnt>?");
            cutSiteGreaterThanPS.setInt(1, 3);// having > 3 tags
            rs = cutSiteGreaterThanPS.executeQuery();
            if (rs.next()) {
            	System.out.println("Number of cut sites with more than 3 tags: " + rs.getInt("numSites"));
            }           
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }

    
    //@Override
    public void putSNPQualityProfile(Map<Position, Map<String,Double>> tagAnnotatedPositionMap, String taxaSubset) {
        int batchCount=0;
        try {
            putSNPPositionsIfAbsent(tagAnnotatedPositionMap.keySet());
            connection.setAutoCommit(false);
            for (Map.Entry<Position, Map<String,Double>> entry : tagAnnotatedPositionMap.entrySet()) {
                Position p=entry.getKey();
                Map<String,Double> vals=entry.getValue();
                int ind=1;

//                snpQualityInsertPS=connection.prepareStatement(
//                        "INSERT OR IGNORE into snpQuality (snpid, taxasubset ,avgDepth, minorDepthProp, minor2DepthProp, gapDepthProp, " +
//                                "propCovered, propCovered2, taxaCntWithMinorAlleleGE2,minorAlleleFreq, inbredF_DGE2)" +
//                                " values(?,?,?,?,?,?,?,?,?,?,?)");
//                System.out.println("vals = " + vals.size());
//                System.out.println(vals.get("inbredF_DGE2").toString());
                snpQualityInsertPS.setInt(ind++, snpPosToIDMap.get(p));
                snpQualityInsertPS.setString(ind++, taxaSubset);
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("avgDepth",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("minorDepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("minor2DepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("gapDepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("propCovered",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("propCovered2",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("taxaCntWithMinorAlleleGE2",0.0));
//                System.out.println("MAF:"+vals.getOrDefault("minorAlleleFreqGE2",-1.0));
                if(vals.containsKey("minorAlleleFreqGE2")) {
                    snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("minorAlleleFreqGE2",0.0));
                    snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("inbredF_DGE2",null));
                }
                else {
                    snpQualityInsertPS.setDouble(ind++,0.0);
                    snpQualityInsertPS.setDouble(ind++,Double.NaN);
                }
                
                //snpQualityInsertPS.get();
                //System.out.println(posTagInsertPS.toString());
                snpQualityInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                    System.out.println("putSNPQualityProfile next"+batchCount);
                    snpQualityInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            snpQualityInsertPS.executeBatch();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }

    public void putSNPQualityProfile(Map<Position, Map<String,Double>> tagAnnotatedPositionMap, String taxaSubset,int counter) throws SQLException {
        try {
            putSNPPositionsIfAbsent(tagAnnotatedPositionMap.keySet());
            connection.setAutoCommit(false);
            for (Map.Entry<Position, Map<String,Double>> entry : tagAnnotatedPositionMap.entrySet()) {
                Position p=entry.getKey();
                Map<String,Double> vals=entry.getValue();
                int ind=1;

//                snpQualityInsertPS=connection.prepareStatement(
//                        "INSERT OR IGNORE into snpQuality (snpid, taxasubset ,avgDepth, minorDepthProp, minor2DepthProp, gapDepthProp, " +
//                                "propCovered, propCovered2, taxaCntWithMinorAlleleGE2,minorAlleleFreq, inbredF_DGE2)" +
//                                " values(?,?,?,?,?,?,?,?,?,?,?)");
//                System.out.println("vals = " + vals.size());
//                System.out.println(vals.get("inbredF_DGE2").toString());
                snpQualityInsertPS.setInt(ind++, snpPosToIDMap.get(p));
                snpQualityInsertPS.setString(ind++, taxaSubset);
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("avgDepth",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("minorDepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("minor2DepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("gapDepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("propCovered",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("propCovered2",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("taxaCntWithMinorAlleleGE2",0.0));
//                System.out.println("MAF:"+vals.getOrDefault("minorAlleleFreqGE2",-1.0));
                if(vals.containsKey("minorAlleleFreqGE2")) {
                    snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("minorAlleleFreqGE2",0.0));
                    snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("inbredF_DGE2",null));
                }
                else {
                    snpQualityInsertPS.setDouble(ind++,0.0);
                    snpQualityInsertPS.setDouble(ind++,Double.NaN);
                }
                if(counter==-1) {
                    snpQualityInsertPS.executeBatch();
                    connection.setAutoCommit(true);    
                }
                else {
                    //snpQualityInsertPS.get();
                    //System.out.println(posTagInsertPS.toString());
                    snpQualityInsertPS.addBatch();
                    if(counter%10000==0) {
                        //System.out.println("putSNPQualityProfile next"+batchCount);
                        snpQualityInsertPS.executeBatch();
                    }
                }
            }
        } catch (SQLException e) {           
            e.printStackTrace();
            throw e;
        }

    }
    private int getMappingApproachID(Position p) throws SQLException{
        String mapApp=p.getAnnotation().getTextAnnotation("mappingapproach")[0];
        if(mapApp==null) return mappingApproachToIDMap.get("unknown");
        Integer val=mappingApproachToIDMap.get(mapApp);
        if(val==null) {
            connection.createStatement().executeUpdate("insert into mappingApproach (approach, software, protocols) " +
                    "values('"+mapApp+"','unknown','unknown')");
            loadMappingApproachHash();
            return mappingApproachToIDMap.get(mapApp);
        } else return val;
    }

    @Override
    public void setTagAlignmentBest(Tag tag, Position position, boolean isBest) {

    }

    @Override
    public boolean putTagAlleles(Multimap<Tag, Allele> tagAlleleMap) {
        int batchCount=0;
        try {
            PreparedStatement alleleTagInsertPS=connection.prepareStatement(
                    "INSERT OR IGNORE into tagallele (alleleid, tagid) values(?,?)");
            putAllTag(tagAlleleMap.keySet());
            loadSNPPositionHash(false);
            putSNPPositionsIfAbsent(tagAlleleMap.values().stream()
                    .map(a -> a.position())
                    .distinct()
                    .collect(Collectors.toSet()));
            putAlleleIfAbsent(tagAlleleMap.values().stream()
                    .distinct()
                    .collect(Collectors.toSet()));
            connection.setAutoCommit(false);
            for (Map.Entry<Tag, Allele> tagAlleleEntry : tagAlleleMap.entries()) {
                int ind=1;
                alleleTagInsertPS.setInt(ind++, alleleToIDMap.get(tagAlleleEntry.getValue()));
                alleleTagInsertPS.setInt(ind++, tagTagIDMap.get(tagAlleleEntry.getKey()));
                alleleTagInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                    System.out.println("alleleTagInsertPS next"+batchCount);
                    alleleTagInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            alleleTagInsertPS.executeBatch();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    @Override
    public boolean putTagAlignmentApproach(String tagAlignmentName, String protocol) {
        return false;
    }

    @Override
    public TaxaDistribution getTaxaDistribution(Tag tag) {
        int tagid=tagTagIDMap.get(tag);
        try {
            tagTaxaDistPS.setInt(1,tagid);
            ResultSet rs=tagTaxaDistPS.executeQuery();
            return TaxaDistBuilder.create(rs.getBytes(1));
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return null;
    }

    @Override
    public Set<Allele> getAlleles(Tag tag) {if(alleleToIDMap==null) loadAlleleHash();
        ImmutableSet.Builder<Allele> alleleBuilder=new ImmutableSet.Builder<>();
        try{
            tagAlleleWhereTagPS.setInt(1,tagTagIDMap.get(tag));
            ResultSet rs=tagAlleleWhereTagPS.executeQuery();
            while(rs.next()) {
                alleleBuilder.add(alleleToIDMap.inverse().get(rs.getInt("alleleid")));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return alleleBuilder.build();
    }

    @Override
    public Multimap<Tag, Allele> getAlleleMap() {
        //if slow consider caching the hash and equal codes
        if(alleleToIDMap==null) loadAlleleHash();
        ImmutableMultimap.Builder<Tag, Allele> tagAlleleBuilder=new ImmutableMultimap.Builder<>();
        try{
            ResultSet rs=connection.createStatement().executeQuery("select * from tagallele");;
            while(rs.next()) {
                tagAlleleBuilder.put(tagTagIDMap.inverse().get(rs.getInt("tagid")), alleleToIDMap.inverse().get(rs.getInt("alleleid")));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return tagAlleleBuilder.build();
    }

    @Override
    public Set<Tag> getTagsForSNPPosition(Position position) {
        ImmutableSet.Builder<Tag> tagBuilder=new ImmutableSet.Builder<>();
        try{
            tagidWhereSNPidPS.setInt(1,snpPosToIDMap.get(position));
            ResultSet rs=tagAlleleWhereTagPS.executeQuery();
            while(rs.next()) {
                tagBuilder.add(tagTagIDMap.inverse().get(rs.getInt("tagid")));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return tagBuilder.build();
    }

    @Override
    public Set<Tag> getTagsForAllele(Position position, byte allele) {
        return getTagsForAllele(new SimpleAllele(allele,position));
    }

    @Override
    public Set<Tag> getTagsForAllele(Allele allele) {
        return null;
    }

    public Multimap<Allele,TaxaDistribution> getAllelesTaxaDistForSNP(Position position) {
        ImmutableMultimap.Builder<Allele,TaxaDistribution> atdBuilder=ImmutableMultimap.builder();
        try{
            alleleTaxaDistForSnpidPS.setInt(1, snpPosToIDMap.get(position));
            ResultSet rs= alleleTaxaDistForSnpidPS.executeQuery();
            while(rs.next()) {
                Allele allele=new SimpleAllele((byte)rs.getInt("allelecall"),position);
                atdBuilder.put(allele,TaxaDistBuilder.create(rs.getBytes("depthsRLE")));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return atdBuilder.build();
    }
    
    @Override
    public Multimap<Allele, Map<Tag, TaxaDistribution>> getAllelesTagTaxaDistForSNP(Position position) {
    	ImmutableMultimap.Builder<Allele,Map<Tag,TaxaDistribution>> aTTdBuilder = ImmutableMultimap.builder();
        
        try{
            if(snpPosToIDMap==null) {
                loadSNPPositionHash(false);
            }      	
            alleleTaxaDistForSnpidPS.setInt(1, snpPosToIDMap.get(position));
            ResultSet rs= alleleTaxaDistForSnpidPS.executeQuery();
            while(rs.next()) {
                Allele allele=new SimpleAllele((byte)rs.getInt("allelecall"),position);
                Tag snpTag = tagTagIDMap.inverse().get(rs.getInt("tagid"));
                TaxaDistribution snpTD = TaxaDistBuilder.create(rs.getBytes("depthsRLE"));
                ImmutableMap.Builder<Tag,TaxaDistribution> tagTDBuilder = ImmutableMap.builder();
                tagTDBuilder.put(snpTag, snpTD);
                aTTdBuilder.put(allele, tagTDBuilder.build());
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return aTTdBuilder.build();
    }

    public Stream<ImmutableMultimap<Allele,TaxaDistribution>> getAllAllelesTaxaDistForSNP() {
        if(snpPosToIDMap==null) {
            loadSNPPositionHash(false);
        }
        Stream<ImmutableMultimap<Allele,TaxaDistribution>> stream = SQL.stream(connection, "select a.*, td.* from allele a, tagallele ta, tagtaxadistribution td\n" +
                "where a.alleleid=ta.alleleid and ta.tagid=td.tagid order by a.snpid")
                .map(entry -> {
                    ImmutableMultimap.Builder<Allele,TaxaDistribution> atdBuilder = ImmutableMultimap.builder();
                    Position pos = snpPosToIDMap.inverse().get(entry.asInt("snpid"));
                    Allele allele=new SimpleAllele((byte)entry.asInt("allelecall"),pos);
                    byte[] byteArray = (byte[])entry.val("depthsRLE").get();
                    
                    atdBuilder.put(allele,TaxaDistBuilder.create(byteArray));
                    return atdBuilder.build();
                });
        return stream;
    }
    public Stream<Map.Entry<Allele,TaxaDistribution>> getAllAllelesTaxaDistForSNPEntries() {
        if(snpPosToIDMap==null) {
            loadSNPPositionHash(false);
        }
        Stream<Map.Entry<Allele,TaxaDistribution>> stream = SQL.stream(connection, "select a.snpid, a.allelecall, td.depthsRLE from allele a, tagallele ta, tagtaxadistribution td\n" +
                "where td.tagid = ta.tagid AND a.alleleid = ta.alleleid order by a.snpid")
                .map(entry -> {
                    Position pos = snpPosToIDMap.inverse().get(entry.asInt("snpid"));
                    Allele allele=new SimpleAllele((byte)entry.asInt("allelecall"),pos);
                    byte[] byteArray = (byte[])entry.val("depthsRLE").get();
                    return new AbstractMap.SimpleEntry(allele,TaxaDistBuilder.create(byteArray));
                });
        return stream;
    }

    @Override
    public Phenotype getAllCountsForTagTissue(Tag tag, String tissue) {
        Integer tagID=tagTagIDMap.get(tag);
        Integer tissueID=tissueTissueIDMap.get(tissue);
        if (tagID == null || tissueID == null) {
            System.out.println("getAllCountsForTagTissue: ERROR, either tissueID or tagID not found in DB");
            return null;
        }
        List<Taxon> taxaList=new ArrayList<>();
        FloatArrayList countList=new FloatArrayList();
        SQL.stream(connection,
                "select t.taxonid, t.readCount from tagTaxaTissueDist t where t.tagid="+tagID+" AND t.tissueid="+tissueID+";")
                .forEach(entry -> {
                    taxaList.add(myTaxaList.get(entry.asInt("taxonid")));
                    countList.add(entry.asInt("readCount"));
                });
        NumericAttribute trait = new NumericAttribute(tissue+":"+tag.sequence()+":count", countList.elements(), new OpenBitSet(countList.size()));
        return new PhenotypeBuilder().fromAttributeList(
                Arrays.asList(new TaxaAttribute(taxaList),trait),
                Arrays.asList(Phenotype.ATTRIBUTE_TYPE.taxa, Phenotype.ATTRIBUTE_TYPE.data)).build().get(0);
    }

    @Override
    public TableReport getAllCountsForTaxonTissue(Taxon taxon, String tissue) {
        // User may query a tissue type or taxon that is not present in the DB
        Integer taxonID=myTaxaList.indexOf(taxon.getName());
        Integer tissueID=tissueTissueIDMap.get(tissue);
        if (tissueID == null || taxonID == null) {
            System.out.println("getAllCountsForTaxonTissue: ERROR, either tissueID or taxonID not found in DB");
            return null;
        }
        String[] headers={"Sequence","TagName","ReadCount"};
        TableReportBuilder reportBuilder=TableReportBuilder.getInstance("CountsFor:"+taxon.getName()+":"+tissue,headers);
        SQL.stream(connection,
                "select t.tagid, t.readCount, tag.tagName from tagTaxaTissueDist t, tag where t.taxonid="+taxonID+" AND t.tissueid="+tissueID
                        +" AND tag.tagid=t.tagid;")
                .forEach(entry -> {
                    reportBuilder.addElements(tagTagIDMap.inverse().get(entry.asInt("tagid")).sequence(),
                            entry.val("tagName").orElse("").toString(), entry.asInt("readCount"));
                });
       return reportBuilder.build();
    }

    @Override
    public Phenotype getAllCountsForTissue(String tissue) {
        Integer tissueID=tissueTissueIDMap.get(tissue);
        if (tissueID == null) {
            System.out.println("getAllCountsForTissue: ERROR, tissue " + tissue + " not found in DB");
            return null;
        }

        final AtomicInteger counter=new AtomicInteger();
        Map<Integer,Integer> taxaIDtoNumAttIndex=SQL.stream(connection, ("select DISTINCT taxonid from tagTaxaTissueDist where tissueid=" + tissueID + " ORDER BY taxonid;"))
                .map(entry -> new Tuple<>(entry.asInt("taxonid"), counter.getAndIncrement()))
                .collect(Collectors.toMap(Tuple::getX, Tuple::getY, (a, b) -> a, TreeMap::new));
        TaxaList taxaList=taxaIDtoNumAttIndex.keySet().stream().map(myTaxaList::get).collect(TaxaList.collect());

        counter.set(0);
        Map<Integer,Integer> tagIDtoNumAttIndex=SQL.stream(connection, ("select DISTINCT tagid from tagTaxaTissueDist where tissueid=" + tissueID + " ORDER BY tagid;"))
                .map(entry -> new Tuple<>(entry.asInt("tagid"),counter.getAndIncrement()))
                .collect(Collectors.toMap(Tuple::getX, Tuple::getY, (a,b) -> a, TreeMap::new));

        float[][] data=new float[tagIDtoNumAttIndex.size()][taxaIDtoNumAttIndex.size()];

        SQL.stream(connection,
                "select t.tagid, t.taxonid, t.readCount from tagTaxaTissueDist t where t.tissueid="+tissueID+";")
                .forEach(entry -> {
                    data[tagIDtoNumAttIndex.get(entry.asInt("tagid"))][taxaIDtoNumAttIndex.get(entry.asInt("taxonid"))]=entry.asInt("readCount");
                });
        List<PhenotypeAttribute> phenotypeAttributes=new ArrayList<>(tagIDtoNumAttIndex.size()+1);
        List<Phenotype.ATTRIBUTE_TYPE> attributeTypeList=new ArrayList<>(tagIDtoNumAttIndex.size()+1);
        phenotypeAttributes.add(new TaxaAttribute(taxaList));
        attributeTypeList.add(Phenotype.ATTRIBUTE_TYPE.taxa);
        tagIDtoNumAttIndex.entrySet().stream().forEachOrdered(entry -> {
            Tag tag = tagTagIDMap.inverse().get(entry.getKey());
            NumericAttribute trait = new NumericAttribute(tissue + ":" + tag.name() + ":count", data[entry.getValue()], new OpenBitSet(data[entry.getValue()].length));
            phenotypeAttributes.add(trait);
            attributeTypeList.add(Phenotype.ATTRIBUTE_TYPE.data);
        });

        return new PhenotypeBuilder().fromAttributeList(phenotypeAttributes,attributeTypeList).build().get(0);
    }

    @Override
    public Set<Tag> getTags() {
        return tagTagIDMap.keySet();
    }

	@Override
    public PositionList getSNPPositions() {
        if(snpPosToIDMap==null) loadSNPPositionHash(false);
        return new PositionListBuilder().addAll(snpPosToIDMap.keySet()).build();
    }

	@Override
    public PositionList getSNPPositions(int minSupportValue) {
		return null;
	}
	
	@Override
    public PositionList getSNPPositions(double minQualityScore) {
		// Add all positions whose quality score equals or exceeds 
		// caller's minQualityScore
		if (minQualityScore == 0) return getSNPPositions();
		if(snpPosToIDMap==null) loadSNPPositionHash(false);
		PositionListBuilder plb = new PositionListBuilder();
		snpPosToIDMap.keySet().stream()
		.forEach(pos -> {
			double[] qs = pos.getAnnotation().getQuantAnnotation("QualityScore");
			if (qs != null && qs.length > 0) {
				if (qs[0] >=minQualityScore) {
					plb.add(pos);
				}
			}			
		});
		return plb.build();
	}
	
    @Override
    public Set<Tag> getTagsForTaxon(Taxon taxon) {
        ImmutableSet.Builder<Tag> tagBuilder=new ImmutableSet.Builder<>();
        int taxonIndex=myTaxaList.indexOf(taxon);
        try {
            ResultSet rs=connection.createStatement().executeQuery("select * from tagtaxadistribution");
            while(rs.next()) {
                if(TaxaDistBuilder.create(rs.getBytes("depthsRLE")).depths()[taxonIndex]>0) {
                    tagBuilder.add(tagTagIDMap.inverse().get(rs.getInt("tagid")));
                }
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return tagBuilder.build();
    }

    @Override
    public Map<Tag, Integer> getTagDepth(Taxon taxon, Position position) {
        return null;
    }

    @Override
    public Map<Tag, Integer> getTagsWithDepth(int minimumDepth) {
        ImmutableMap.Builder<Tag, Integer> tagBuilder=new ImmutableMap.Builder<>();
        try {
            ResultSet rs=connection.createStatement().executeQuery(
                    "select tagid, totalDepth from tagtaxadistribution where totalDepth >= "+minimumDepth);
            while(rs.next()) {
                tagBuilder.put(tagTagIDMap.inverse().get(rs.getInt(1)),rs.getInt(2));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return tagBuilder.build();
    }

    @Override
    public PositionList getTagCutPositions(boolean onlyBest) {
        if(cutPosToIDMap == null) loadCutPositionHash();
        PositionListBuilder plb=new PositionListBuilder();
        cutPosToIDMap.keySet().stream()
                //.filter(p -> p.isAnnotatedWithValue("isbest","true"))
                .forEach(p -> plb.add(p));  //todo only best not implemented here
        plb.sortPositions();
        return plb.build();
    }

    @Override
    public PositionList getTagCutPositions(Chromosome chromosome, int firstPosition, int lastPosition, boolean onlyBest) {
        PositionListBuilder plb=new PositionListBuilder();
        plb.addAll(getPositionSubMap(chromosome,firstPosition,lastPosition).keySet());  //todo only best not implemented here
        return plb.build();
    }

 
    /**
     * Get the cut position associated with each tag in a set.  Return a map of Tag/Position
     * from which the cut position/strand will be pulled.
     * This is used in debugging with the SNPCutPosTagVerificationPlugin
     * @return map of Tag/Position 
     */
    @Override
    public Map<Tag, Position> getTagCutPosition(Set<Tag> tagSet){
    	ImmutableMap.Builder<Tag,Position> tagPosBuilder=new ImmutableMap.Builder<>();
    	for (Tag tag : tagSet) {
    		int tagID=tagTagIDMap.get(tag);       
    		try {
    			ResultSet rs=connection.createStatement()
    					.executeQuery("select cp.*, tcp.* from cutposition cp, tag t, tagCutPosition tcp " +
    							"where tcp.tagid=t.tagid and tcp.positionid=cp.positionid and t.tagid= " + tagID);
    			while(rs.next()) { // create the position, add to map   				
    				Position cutPos=new GeneralPosition
    						.Builder(new Chromosome(rs.getString("chromosome")),rs.getInt("position"))
    				.strand(rs.getByte("strand"))
    				.addAnno("forward", rs.getBoolean("forward")?"true":"false")
    				.build();
    				tagPosBuilder.put(tag,cutPos);
    			}               
    		} catch (SQLException e) {
    			e.printStackTrace();
    			return null;
    		}
    	}
    	return tagPosBuilder.build();
    }
    
    private Map<Position,Integer> getPositionSubMap(Chromosome chromosome, int firstPosition, int lastPosition) {
        if(cutPosToIDMap==null) loadCutPositionHash();
        Position startPos=new GeneralPosition.Builder(chromosome,firstPosition).build();
        if(lastPosition<0) lastPosition=Integer.MAX_VALUE;
        Position lastPos=new GeneralPosition.Builder(chromosome,lastPosition).build();
        return cutPosToIDMap.subMap(startPos,lastPos);
    }

    @Override
    public Map<String, String> getTagAlignmentApproaches() {
        ImmutableMap.Builder<String,String> appBuilder=new ImmutableMap.Builder<>();
        try {
            ResultSet rs=connection.createStatement().executeQuery("select * from mappingApproach");
            while(rs.next()) {
                appBuilder.put(rs.getString("approach"), rs.getString("software") + ":" + rs.getString("approach"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return appBuilder.build();
    }

//    @Override
//    public Map<Position, Map<Tag, TaxaDistribution>> getCutPositionTagTaxaMapX(Chromosome chromosome, int firstPosition, int lastPosition) {
//        //consider doing this all by SQL if performance suffers
//        PositionList pl=getTagCutPositions(chromosome,firstPosition,lastPosition,true);
//        ImmutableMap.Builder<Position, Map<Tag, TaxaDistribution>> positionMapBuilder=new ImmutableMap.Builder<>();
//        pl.stream().forEach(p -> positionMapBuilder.put(p,getTagsTaxaMap(p)));
//        //this is slow as each position is a separate transaction
//        return positionMapBuilder.build();
//    }

    //TODO need to add the forward direction to the resulting map somehow.  Perhaps  Map<Position, Map<Tag, Tuple<Boolean,TaxaDistribution>>>
    //alternatively there could be a tag alignment object.

    @Override
    public Map<Position, Map<Tag, Tuple<Boolean,TaxaDistribution>>> getCutPositionTagTaxaMap(Chromosome chromosome, int firstPosition, int lastPosition) {
        // With the change to process Chromosomes as strings instead of int,
        // this query throws an "SQL error or missing database (unrecognized token: "7D")"
        // error when the chrom name contains non-digit characters.  The chromosome field in the 
        // cutposition table is declared as "TEXT" but without the single quotes, SQL takes the
        // input value as an integer.  Adding single quotes around the string fixes the problem.
        String sqlQuery="select p.positionid, forward, chromosome, position, strand, t.tagid, depthsRLE  " +
                "from tag t, cutposition p, tagCutPosition tc, tagtaxadistribution ttd " +
                "where p.positionid=tc.positionid and tc.tagid=t.tagid and t.tagid=ttd.tagid " +
                "and chromosome='"+chromosome.toString()+"'" +//" and position>"+firstPosition+" " + //todo position would need to be index to make fast
                " order by position";
        Map<Position, Map<Tag, Tuple<Boolean,TaxaDistribution>>> positionTagTaxaMap=new HashMap<>();
        Map<Integer,Position> tempPositionMap=new HashMap<>();  //reverse the map
        getPositionSubMap(chromosome,firstPosition,lastPosition).entrySet().stream()
                .forEach(entry -> tempPositionMap.put(entry.getValue(),entry.getKey()));
        // Adding totalTagsAdded so this method can be used to verify the counts when
        // the new method getCutPosForSTrandTagTaxaMap() is called for forward and
        // reverse strands.  
        int totalTagsAdded = 0; 
        try{
            ResultSet rs=connection.createStatement().executeQuery(sqlQuery);
            while(rs.next()) {
                Position position=tempPositionMap.get(rs.getInt("positionid"));
                Tag tag=tagTagIDMap.inverse().get(rs.getInt("tagid"));
                TaxaDistribution taxaDistribution=TaxaDistBuilder.create(rs.getBytes("depthsRLE"));
                Boolean forwardAlignDirection=rs.getBoolean("forward");
                Map<Tag, Tuple<Boolean,TaxaDistribution>> tagTaxaMap=positionTagTaxaMap.computeIfAbsent(position, k -> new HashMap<>());
                tagTaxaMap.put(tag,new Tuple<>(forwardAlignDirection,taxaDistribution));
                totalTagsAdded++;
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
       
        System.out.println("positionTagTaxaMap = " + positionTagTaxaMap.size() +
                " totalTagsAdded: " + totalTagsAdded);
        return positionTagTaxaMap;
    }

    @Override
    public Map<Position, Map<Tag, TaxaDistribution>> getCutPosForStrandTagTaxaMap(Chromosome chromosome, int firstPosition, int lastPosition, boolean direction) {
        // With the change to process Chromosomes as strings instead of int,
        // this query throws an "SQL error or missing database (unrecognized token: "7D")"
        // error when the chrom name contains non-digit characters.  The chromosome field in the 
        // cutposition table is declared as "TEXT" but without the single quotes, SQL takes the
        // input value as an integer.  Adding single quotes around the string fixes the problem.
        String sqlQuery="select p.positionid, forward, chromosome, position, strand, t.tagid, depthsRLE  " +
                "from tag t, cutposition p, tagCutPosition tc, tagtaxadistribution ttd " +
                "where p.positionid=tc.positionid and tc.tagid=t.tagid and t.tagid=ttd.tagid " +
                "and chromosome='"+chromosome.toString()+"'" +//" and position>"+firstPosition+" " + //todo position would need to be index to make fast
                " order by position";

        Map<Position, Map<Tag, TaxaDistribution>> positionTagTaxaMap=new HashMap<>();
        Map<Integer,Position> tempPositionMap=new HashMap<>();  //reverse the map
        getPositionSubMap(chromosome,firstPosition,lastPosition).entrySet().stream()
                .forEach(entry -> tempPositionMap.put(entry.getValue(),entry.getKey()));
        int totalTagsAdded = 0;
        try{
            ResultSet rs=connection.createStatement().executeQuery(sqlQuery);
            while(rs.next()) {              
                Boolean forwardAlignDirection=rs.getBoolean("forward");
                if (forwardAlignDirection == direction) {                   
                    // tags from forward and reverse strands will be aligned separately
                    Position position=tempPositionMap.get(rs.getInt("positionid"));
                    Tag tag=tagTagIDMap.inverse().get(rs.getInt("tagid"));
                    TaxaDistribution taxaDistribution=TaxaDistBuilder.create(rs.getBytes("depthsRLE"));
                    Map<Tag, TaxaDistribution> tagTaxaMap=positionTagTaxaMap.computeIfAbsent(position, k -> new HashMap<>());
                    tagTaxaMap.put(tag,taxaDistribution);
                    totalTagsAdded++;
                }               
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        // The positionsTagTaxaMap size is based on number of positions, not number of tags.
        // Tags added should add up to number of forward tags plus number of reverse tags.
        System.out.println("positionTagTaxaMap size, forwardStrand, " + direction + " size:" + positionTagTaxaMap.size()
          +  " totalTagsAdded:" + totalTagsAdded);
        return positionTagTaxaMap;
    }

    @Override
    public Map<Tag, TaxaDistribution> getTagsTaxaMap(Position cutPosition) {
        ImmutableMap.Builder<Tag, TaxaDistribution> tagTaxaDistributionBuilder=new ImmutableMap.Builder<>();
        if(cutPosToIDMap==null) loadCutPositionHash();
        try{
            taxaDistWhereCutPositionIDPS.setInt(1,cutPosToIDMap.get(cutPosition));
            ResultSet rs=taxaDistWhereCutPositionIDPS.executeQuery();
            while(rs.next()) {
                tagTaxaDistributionBuilder.put(tagTagIDMap.inverse().get(rs.getInt("tagid")), TaxaDistBuilder.create(rs.getBytes("depthsRLE")));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return tagTaxaDistributionBuilder.build();
    }

    private void putCutPositionsIfAbsent(Collection<Position> positions) {
        try {
        int batchCount=0;
        if(cutPosToIDMap==null) loadCutPositionHash();
        connection.setAutoCommit(false);
        PreparedStatement posInsertPS=connection.prepareStatement(
                "INSERT OR IGNORE into cutposition (chromosome, position, strand) values(?,?,?)");
        for (Position p : positions) {
            if(cutPosToIDMap.containsKey(p)) continue;
            posInsertPS.setString(1, p.getChromosome().toString());
            posInsertPS.setInt(2, p.getPosition());
            posInsertPS.setByte(3, p.getStrand());
            posInsertPS.addBatch();
            batchCount++;
            if(batchCount>10000) {
                System.out.println("putCutPositionsIfAbsent next"+batchCount);
                posInsertPS.executeBatch();
                batchCount=0;
            }
        }
        posInsertPS.executeBatch();
        if(batchCount>0) loadCutPositionHash();
        connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }

    private void putSNPPositionsIfAbsent(Collection<Position> positions) {
        try {
            int batchCount=0;
            if(snpPosToIDMap==null) loadSNPPositionHash(false);
            connection.setAutoCommit(false);
            PreparedStatement snpPosInsertPS=connection.prepareStatement(
                    "INSERT OR IGNORE into snpposition (chromosome, position, strand,qualityScore,refAllele) values(?,?,?,?,?)");
            for (Position p : positions) {
                if(snpPosToIDMap.containsKey(p)) continue;
                snpPosInsertPS.setString(1, p.getChromosome().toString());
                snpPosInsertPS.setInt(2, p.getPosition());
                snpPosInsertPS.setByte(3, p.getStrand());
                double[] qsList = p.getAnnotation().getQuantAnnotation("QualityScore");               
                if (qsList != null & qsList.length > 0) {
                	snpPosInsertPS.setFloat(4, (float)qsList[0]);
                } else {
                	snpPosInsertPS.setFloat(4, (float)0.0);
                } 
                // Position is now annotated with reference, add it
                // from the position's stored Reference - should be created
                // in Discovery's referencePositions() method
                snpPosInsertPS.setByte(5,p.getAllele(WHICH_ALLELE.Reference));
 
                snpPosInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                    System.out.println("putSNPPositionsIfAbsent next"+batchCount);
                    snpPosInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            snpPosInsertPS.executeBatch();
            if(batchCount>0) loadSNPPositionHash(false);
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void putAlleleIfAbsent(Collection<Allele> alleles) throws IllegalStateException {
        try {
            int batchCount=0;
            PreparedStatement alleleInsertPS=connection.prepareStatement(
                    "INSERT OR IGNORE into allele (snpid, allelecall, qualityscore) values(?,?,?)");
            connection.setAutoCommit(false);
            for (Allele allele : alleles) {
                Integer snpID=snpPosToIDMap.get(allele.position());
                if(snpID==null) throw new IllegalStateException("SNP position missing for allele");
                int index=1;
                alleleInsertPS.setInt(index++, snpID);
                alleleInsertPS.setByte(index++, allele.allele());
                alleleInsertPS.setByte(index++, (byte) 0);  //todo set quality scores once annotation pattern is set
                alleleInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                    System.out.println("putAlleleIfAbsent next"+batchCount);
                    alleleInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            alleleInsertPS.executeBatch();
            if(batchCount>0) loadAlleleHash();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }
    
    // Return SNP positions for specified list of chromosomes
    @Override
    public  PositionList getSNPPositionsForChromosomes(Integer startChr, Integer endChr) {
    	PositionListBuilder plb = new PositionListBuilder();
    	// Verify good chromsome values
    	if (startChr < 1 ||
    	    endChr < 1 ||
    	    startChr > endChr ) {
    		System.err.printf("getSNPPOsitionsForChromosomes:  bad Chromosome values: startChr %d, endChr %d\n",
    				startChr, endChr);
    		return null;
    	}
        try{
        	for (int chrom = startChr; chrom <= endChr; chrom++ ){
        		snpPositionsForChromosomePS.setString(1, Integer.toString(chrom));
                ResultSet rs=snpPositionsForChromosomePS.executeQuery();
                while(rs.next()) {
                	Chromosome chr= new Chromosome(Integer.toString(chrom));
                	byte refAllele = (byte)rs.getInt("refAllele");
                	Position position = new GeneralPosition.Builder(chr,rs.getInt("position"))
                	.addAnno("QualityScore",rs.getFloat("qualityScore"))
                	.allele(WHICH_ALLELE.Reference,refAllele).build();
                    plb.add(position);
                }
        	} 
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return plb.build();
    }

    // Store SNP quality score for positions in snpposition table
	@Override
	public void putSNPPositionQS(PositionList qsPL) {
		int batchCount=0;
		try {
			connection.setAutoCommit(false);
			PreparedStatement qsUpdatePS = connection.prepareStatement("update snpposition set qualityScore = ? where chromosome = ? and position = ?");
			// if row to update doesn't exist, update skips and continues.  A 0 return
			// for any entry in qsUpdatePS.executeBatch() will mean update wasn't performed.
			for (Position qsPos : qsPL) {
	            double qscore = 0.0;
				double[] qs = qsPos.getAnnotation().getQuantAnnotation("QualityScore");
				if (qs != null && qs.length > 0) {
					qscore = qs[0];	
				}			
	            qsUpdatePS.setFloat(1, (float)qscore);
	            qsUpdatePS.setString(2,qsPos.getChromosome().getName());
	            qsUpdatePS.setInt(3,qsPos.getPosition());
	            qsUpdatePS.addBatch();
	            batchCount++;
				if(batchCount>100000) {
					System.out.println("updateSNPPosition next "+batchCount);
					qsUpdatePS.executeBatch();
					batchCount=0;
				}
			}				
			//int[] lastBatch = qsUpdatePS.executeBatch(); // use to see results
			qsUpdatePS.executeBatch();
			connection.setAutoCommit(true);
		} catch (SQLException e) {
			System.out.println("Error executing UPDATE statement for TagDataSQLite:putSNPQualityPositions");
			e.printStackTrace();
		}
		loadSNPPositionHash(true); // reload this map to obtain the new position quality scores.
	}

	// THis is intended to be used for junit verification, where there the query
	// numbers will be small.
	//  The method will need to be modified with batching if it is to be used in a more
	// query intensive environment.
	@Override
	public ListMultimap<String, Tuple<Integer, Float>> getSNPPositionQS(HashMultimap<String, Integer> myMap) {
	    ImmutableListMultimap.Builder<String, Tuple<Integer, Float>> qsMap = new ImmutableListMultimap.Builder<String, Tuple<Integer, Float>>()
	                .orderKeysBy(Ordering.natural()).orderKeysBy(Ordering.natural());
	    int batchCount = 0;
		try {
			connection.setAutoCommit(false);
			PreparedStatement dbTestPS = connection.prepareStatement("select qualityScore from snpposition where chromosome = ? and position = ?");
			for (Map.Entry<String, Integer> entry : myMap.entries()) {
				String chrom = entry.getKey();
				Integer posQS = entry.getValue();
				dbTestPS.setString(1, chrom);
				dbTestPS.setInt(2, posQS);
				//dbTestPS.addBatch();
				//batchCount++;
				ResultSet rs=dbTestPS.executeQuery();
	            while(rs.next()) {
	                Float myQual = rs.getFloat("qualityScore");
	                Tuple<Integer,Float> myTuple = new Tuple<Integer,Float>(posQS,myQual);
	                qsMap.put(chrom, myTuple);
	            }
			}
		} catch (Exception exc) {
			exc.printStackTrace();
			return null;
		}
		return qsMap.build();
	}

    @Override
    public List<Chromosome> getChromosomesFromCutPositions() {
        // Get a list of chromosomes from the cutPosition table
        List<Chromosome> chromList = new ArrayList<>();
        try {
            ResultSet rs = connection.createStatement().executeQuery("select DISTINCT chromosome from cutPosition");
            while(rs.next()) {
                Chromosome chrom=new Chromosome(rs.getString("chromosome"));
                chromList.add(chrom);
            }
        } catch (SQLException exc) {
            exc.printStackTrace();
            return null;
        }
        return chromList;
    }

    @Override
    public void clearTagTaxaDistributionData() {
        // Clear tags and tagtaxadist tables populated from GBSSeqToTagDBPlugin
        try {
            boolean rs = connection.createStatement().execute("delete FROM tagtaxadistribution"); 
            rs = connection.createStatement().execute("delete FROM tag");
            tagTagIDMap = null;
            loadTagHash();
        } catch (SQLException exc) {
            System.out.println("ERROR - problem deleting tagtaxadistribution data");
            exc.printStackTrace();
        }              
    }
    
    @Override
    public void clearAlignmentData() {
        // Clear tables populated from SAMToGBSdbPluging
        try {           
            boolean rs = connection.createStatement().execute("delete FROM tagCutPosition");
            rs = connection.createStatement().execute("delete FROM cutPosition");
            rs = connection.createStatement().execute("delete FROM mappingApproach");
            cutPosToIDMap = null;
            mappingApproachToIDMap = null;
            loadMappingApproachHash(); // this adds "unknown" to mappingApproachToIDMap
        } catch (SQLException exc) {
            System.out.println("ERROR - problem deleting alignment data");
            exc.printStackTrace();
        }              
    }
    
    @Override
    public void clearDiscoveryData() {
        // Clear all entries from tables populated from the DiscoverySNPCallerPluginV2 
        // The "delete" removes data, but keeps the table size.  The "vacuum" command
        // is NOT called as it rebuilds the entire data base which can be quite time intensive.
        // The rows needed for this table will be needed again in the subsequent run.  Vacuum also creates
        // a new file, so disk space requirements could double while vacuuming.
        try {         
            boolean rs = connection.createStatement().execute("delete FROM tagallele");            
            rs = connection.createStatement().execute("delete FROM snpposition"); 
            rs = connection.createStatement().execute("delete FROM allele");
            alleleToIDMap = null;
            snpPosToIDMap = null;
        } catch (SQLException exc) {
            System.out.println("ERROR - problem deleting discovery data");
            exc.printStackTrace();
        }       
    }

    @Override
    public void clearSNPQualityData() {
        // Clear table populated via SNPQualityProfilerPlugin
        try {
            boolean rs = connection.createStatement().execute("delete FROM snpQuality");         
        } catch (SQLException exc) {
            System.out.println("ERROR - problem deleting snpQuality data");
            exc.printStackTrace();
        }       
    }
    

    // This method is called from GBSSeqToTagDBPlugin when a user wishes to append data.
    // The map created must NOT be a fixed taxaDist map.  We intend to increment the values
    @Override
    public Map<Tag, TaxaDistribution> getAllTagsTaxaMap() {
        ImmutableMap.Builder<Tag,TaxaDistribution> tagTDBuilder = ImmutableMap.builder();
        loadTagHash(); // get updated tag map
        try {
            ResultSet rs = connection.createStatement().executeQuery("select tagid, depthsRLE from tagtaxadistribution");
            while(rs.next()) {
                Tag myTag = tagTagIDMap.inverse().get(rs.getInt("tagid"));
                TaxaDistribution myTD = TaxaDistBuilder.create(rs.getBytes("depthsRLE"));
                TaxaDistribution myTDIncr = TaxaDistBuilder.create(myTD); // convert to incrementable version
                tagTDBuilder.put(myTag, myTDIncr);
            }
            return tagTDBuilder.build();
        } catch (SQLException exc) {
            System.out.println("getAllTaxaMap: caught SQLException attempting to grab taxa Distribution ");
            exc.printStackTrace();
        }
        return tagTDBuilder.build();
    }

    
    @Override
    public boolean putAllTissue(ArrayList<String> tissues) {
        int batchCount=0, totalCount=0;
        try {
            connection.setAutoCommit(false);
            PreparedStatement tissueInsertPS=connection.prepareStatement("insert into tissue (tissue) values(?)");
            for (String tissue : tissues) {
                if(tissueTissueIDMap.containsKey(tissue)) continue;  //it is already in the DB skip
                tissueInsertPS.setString(1, tissue);
                tissueInsertPS.addBatch();
                batchCount++;
                totalCount++;
                if(batchCount>100000) {
                    System.out.println("tissueInsertPS.executeBatch() "+batchCount);
                    tissueInsertPS.executeBatch();
                    //connection.commit();
                    batchCount=0;
                }
            }
            tissueInsertPS.executeBatch();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
            return false;
        }
        if(totalCount>0) loadTissueHash();
        return true;
    }

    @Override
    public Set<String> getAllTissue() {
        return tissueTissueIDMap.keySet();
    }


}

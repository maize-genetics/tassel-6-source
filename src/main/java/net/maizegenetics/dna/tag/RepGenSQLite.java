
/**
 * @author lcj34
 *
 */

package net.maizegenetics.dna.tag;

import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.sqlite.SQLiteConfig;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultimap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Ordering;
import com.google.common.io.CharStreams;

import net.maizegenetics.analysis.gbs.repgen.AlignmentInfo;
import net.maizegenetics.analysis.gbs.repgen.RefTagData;
import net.maizegenetics.analysis.gbs.repgen.TagCorrelationInfo;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.dna.snp.SimpleAllele;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.db.SQL;


public class RepGenSQLite implements RepGenDataWriter, AutoCloseable {
    private Connection connection = null;

    /*These maps contain  objects that are most queried by users.  This is a not the simplest way to do this, which
        would probably be done cleaner just with queries against the databases.  However, there are large performance
        gains in the case of SQLite (or at least with my ability to optimize).

        The logic behind this most of the datasets are relatively small (except tagTagIDMap), and this prevents creation
        of these objects over and over again.
     */
    private BiMap<Tag,Integer> tagTagIDMap; // Tag is AbstractTag.java
    
    // RefTagData is the key to BiMap refTagRefTagIDMAP because:
    // The reference tag sequence can show up in different places.
    // refTag is unique based on sequence, chromosome, seqlen, position, refGenomeID.  
    // Sequence can be duplicated on same chrom and on multiple chroms.
    private BiMap<RefTagData,Integer> reftagReftagIDMap; // Tag is AbstractRefTag.java
    private Map<String,Integer> mappingApproachToIDMap;
    private BiMap<String,Integer> referenceGenomeToIDMap;
    private SortedMap<Position,Integer> physicalMapPositionToIDMap;
    public BiMap<Position,Integer> snpPosToIDMap;
    private BiMap<Allele,Integer> alleleToIDMap;

    private TaxaList myTaxaList;

    PreparedStatement tagTaxaDistPS;
    PreparedStatement allelePairWithTagid1PS;
    PreparedStatement posTagMappingInsertPS;
    PreparedStatement taxaDistWhereTagMappingIDPS;
    PreparedStatement snpPositionsForChromosomePS;
        
    PreparedStatement snpQualityInsertPS;
    
    PreparedStatement allelePairInsertPS;
    PreparedStatement tagtagAlignmentInsertPS;
    PreparedStatement tagReftagAlignmentInsertPS;
    PreparedStatement refRefAlignmentInsertPS;
    PreparedStatement tagAlignForNonRefTagPS;
    PreparedStatement refAlignForRefTagPS;
    PreparedStatement nonReftagAlignmentsForRefTagPS;
    PreparedStatement refTagAlignsForNonRefTagPS;
    PreparedStatement tagTagCorrelationInsertPS;
    PreparedStatement tagCorrelationsForTag1PS;
    PreparedStatement tagCorrelationsForTag2PS;

    public RepGenSQLite(String filename) {
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

            connection = DriverManager.getConnection("jdbc:sqlite:"+filename,config.toProperties());
            connection.setAutoCommit(true);  //This has massive performance effects
            Statement statement = connection.createStatement();
            statement.setQueryTimeout(30);  // set timeout to 30 sec.
            //                System.out.println(schema);
            if(doesDBExist==false) {
                String schema = CharStreams.toString(new InputStreamReader(RepGenSQLite.class.getResourceAsStream("RepGenSchema.sql")));
                statement.executeUpdate(schema);
            }
            initPreparedStatements();
            loadReferenceGenomeHash();
            loadTagHash();
            loadRefTagHash();
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
            posTagMappingInsertPS=connection.prepareStatement(
                    "INSERT OR IGNORE into TagMapping (reftagid, position_id, method_id, bp_error, cm_error)" +
                    " values(?,?,?,?,?)");
            
            tagTaxaDistPS=connection.prepareStatement("select depthsRLE from tagtaxadistribution where tagid=?");
            allelePairInsertPS = connection.prepareStatement(
                    "INSERT into allelepair (tagid_1,tagid_2,qualityscore)" + 
                    " values(?,?,?)");
            allelePairWithTagid1PS=connection.prepareStatement("select * from allelepair where tagid_1=?");
//            taxaDistWhereTagMappingIDPS=connection.prepareStatement(
//                    "select tagtaxadistribution.* from tagMapping, tagtaxadistribution where tagMapping.position_id=? and " +
//                    "tagMapping.tagid=tagtaxadistribution.tagid");
            snpPositionsForChromosomePS=connection.prepareStatement(
                    "select position, qualityScore, refAllele from snpposition where chromosome=?");

            snpQualityInsertPS=connection.prepareStatement(
                    "INSERT into snpQuality (snpid, taxasubset ,avgDepth, minorDepthProp, minor2DepthProp, gapDepthProp, " +
                            "propCovered, propCovered2, taxaCntWithMinorAlleleGE2, minorAlleleFreqGE2, inbredF_DGE2)" +
                    " values(?,?,?,?,?,?,?,?,?,?,?)");
            //                snpQualityInsertPS=connection.prepareStatement(
            //                        "INSERT into snpQuality (snpid, taxasubset)" +
            //                                " values(?,?)");
            
            tagtagAlignmentInsertPS=connection.prepareStatement(
                    "INSERT into tag_tag_Alignments (tag1id, tag2id, score )" +
                    " values(?,?,?)");
            tagReftagAlignmentInsertPS = connection.prepareStatement(
                    "INSERT into tag_reftag_Alignments(tag1id, refTagID, score,ref_align_start_pos,ref_align_strand )" +
                    " values(?,?,?,?,?)");
            refRefAlignmentInsertPS = connection.prepareStatement(
                    "INSERT into reftag_reftag_Alignments(tag1id, tag2id, score )" +
                    " values(?,?,?)");
            tagTagCorrelationInsertPS=connection.prepareStatement(
                    "INSERT into tagCorrelations (tag1id, tag2id, t1t2_pearson, t1t2_spearman, pres_abs_pearson, r2 )" +
                    " values(?,?,?,?,?,?)");
            // because there can be a tagID X in both tag and refTag table, you must
            // specify that this query only wants the values where tag1 is NOT a ref
            // (IE tag1ID comes from the tagTagIDMap)
            // THis gets both the non-ref and ref tag alignments for a particular non-ref tag.
            // NOTE: SQLite has no real boolean field.  The values are stored as 0 and 1
            tagAlignForNonRefTagPS= connection.prepareStatement(
                    "select tag2id,  score " +
                    "from tag_tag_Alignments where tag1id=? and score >= ?");
            // For this query, tag1_isref is true:  the tag1id is a tag from the refTagRefTagIDMap
            refAlignForRefTagPS = connection.prepareStatement(
                    "select tag2id, score " +
                    "from reftag_reftag_Alignments where tag1id=? and score >= ?");
            nonReftagAlignmentsForRefTagPS = connection.prepareStatement(
                    "select tag1id, score, ref_align_start_pos, ref_align_strand " +
                    "from tag_refTag_Alignments where refTagID=? and score >= ?");
            refTagAlignsForNonRefTagPS = connection.prepareStatement(
                    "select refTagID, score, ref_align_start_pos, ref_align_strand " +
                    "from tag_reftag_Alignments where tag1id=? and score >= ?");
            tagCorrelationsForTag1PS = connection.prepareStatement(
                    "select tag2id, t1t2_pearson, t1t2_spearman, pres_abs_pearson, r2 " +
                    "from tagCorrelations where tag1id=?");
            tagCorrelationsForTag2PS = connection.prepareStatement(
                    "select tag1id, t1t2_pearson, t1t2_spearman, pres_abs_pearson, r2 " +
                    "from tagCorrelations where tag2id=?");
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
    
    private void loadRefTagHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from reftag");
            int size=rs.getInt(1);
            System.out.println("size of all tags in reftag table=" + size);
            if(reftagReftagIDMap==null || size/(reftagReftagIDMap.size()+1)>3) reftagReftagIDMap=HashBiMap.create(size);
            rs=connection.createStatement().executeQuery("select * from refTag");
            boolean hasName;
            try{rs.findColumn("tagName");hasName=true;}catch (SQLException e) {hasName=false;}
            while(rs.next()) {
                TagBuilder tagBuilder=TagBuilder.instance(rs.getBytes("sequence"),rs.getShort("seqlen"));
                if(hasName) tagBuilder.name(rs.getString("tagName"));
                Tag refTag = tagBuilder.build();
                String chrom = rs.getString("chromosome");
                int pos = rs.getInt("position");
                int refGenID = rs.getInt("refGenomeID");
                 
                String refGenName = referenceGenomeToIDMap.inverse().get(refGenID);
                RefTagData rtd = new RefTagData(refTag,chrom,pos,refGenName);
                reftagReftagIDMap.putIfAbsent(rtd,rs.getInt("reftagid"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void loadPhysicalMapPositionHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from physicalmapposition");
            int size=rs.getInt(1);
            System.out.println("Before loading new positions, size of all positions in physicalMapPosiiton table="+size);
            if(physicalMapPositionToIDMap==null) {physicalMapPositionToIDMap=new TreeMap<>();}
            else if(size==physicalMapPositionToIDMap.size()) return;
            rs=connection.createStatement().executeQuery("select * from physicalMapPosition");
            while(rs.next()) {
                Position p=new GeneralPosition
                        .Builder(new Chromosome(rs.getString("chromosome")),rs.getInt("physical_position"))
                        .strand(rs.getByte("strand"))
                        .build();
                physicalMapPositionToIDMap.putIfAbsent(p, rs.getInt("posid"));
            }
            rs=connection.createStatement().executeQuery("select count(*) from physicalmapposition");
            size=rs.getInt(1);
            System.out.println("After loading new positions, size of all positions in physicalMapPosiiton table="
            +size + ", size of physicalMapPositionToIDMAP: " + physicalMapPositionToIDMap.size());
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

    private void loadReferenceGenomeHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from referenceGenome");
            int size=rs.getInt(1);
            System.out.println("size of all references in referenceGenome table="+size);
            if(size==0) {
                connection.createStatement().executeUpdate("insert into referenceGenome (refname) " +
                        "values('unknown')");
                size=1;
            }
            referenceGenomeToIDMap=HashBiMap.create(size);
            rs=connection.createStatement().executeQuery("select * from referenceGenome");
            while(rs.next()) {
                referenceGenomeToIDMap.put(rs.getString("refname"), rs.getInt("refid"));
                System.out.println("refence name from referenceGenome: " + rs.getString("refname"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    private void loadMappingApproachHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from mappingApproach");
            int size=rs.getInt(1);
            System.out.println("size of all approaches in mappingApproach table="+size);
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

    // Only when called from RepGenLoadSeqToDBPLugin do we have the number of instances
    // and the quality score.  Other places (from within this file) don't have that data.
    // Populate table based on information sent in.
    @Override
    public boolean putAllTag(Set<Tag>tags,Map<Tag,Tuple<Integer,String>> tagInstanceAverageQS) {
    //public boolean putAllTag(Set<Tag> tags) {
        int batchCount=0, totalCount=0;
        if (tagInstanceAverageQS != null) {
            try {
                connection.setAutoCommit(false);
                PreparedStatement tagInsertPS=
                        connection.prepareStatement("insert into tag (sequence, seqlen,isReference,qualityScore,numTagInstances) values(?,?,?,?,?)");
                for (Map.Entry<Tag, Tuple<Integer,String>> entry : tagInstanceAverageQS.entrySet()) {
                    Tag tag = entry.getKey();
                    if(tagTagIDMap.containsKey(tag)) continue;  //it is already in the DB skip
                    int numInstances = entry.getValue().x;
                    String qscore = entry.getValue().y;
                    tagInsertPS.setBytes(1, tag.seq2BitAsBytes());
                    tagInsertPS.setShort(2, tag.seqLength());
                    tagInsertPS.setBoolean(3, tag.isReference());
                    tagInsertPS.setString(4, qscore);
                    tagInsertPS.setInt(5, numInstances);
                    tagInsertPS.addBatch();
                    batchCount++;
                    totalCount++;
                    if(batchCount>10000) {
                       // System.out.println("tagInsertPS.executeBatch() "+batchCount);
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
        } else {
            try {
                connection.setAutoCommit(false);
                PreparedStatement tagInsertPS=
                        connection.prepareStatement("insert into tag (sequence, seqlen,isReference) values(?,?,?)");
                for (Tag tag: tags) {                
                    if(tagTagIDMap.containsKey(tag)) continue;  //it is already in the DB skip                    
                    tagInsertPS.setBytes(1, tag.seq2BitAsBytes());
                    tagInsertPS.setShort(2, tag.seqLength());
                    tagInsertPS.setBoolean(3, tag.isReference());                   
                    tagInsertPS.addBatch();
                    batchCount++;
                    totalCount++;
                    if(batchCount>100000) {
                       // System.out.println("tagInsertPS.executeBatch() "+batchCount);
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
        }

        if(totalCount>0) {
            System.out.println("RepGenSQLite:putAllTag, totalCount=" + totalCount + ",loadingHash");
            loadTagHash();
        } 
        return true;
    }

    
    // Only when called from RepGenLoadSeqToDBPLugin do we have the number of instances
    // and the quality score.  Other places (from within this file) don't have that data.
    // Populate table based on information sent in.
    @Override
    public boolean putAllRefTag(Multimap<Tag,Position> refTagPositionMap, String refGenome) {
        int batchCount=0, totalCount=0;
        try {
            connection.setAutoCommit(false);
            int refGenomeID = referenceGenomeToIDMap.get(refGenome);
            PreparedStatement refTagInsertPS=
                    connection.prepareStatement("insert into reftag (sequence, seqlen, chromosome,position,refGenomeID) values(?,?,?,?,?)");
            for (Map.Entry<Tag, Position> entry : refTagPositionMap.entries()) {
                Tag tag = entry.getKey();
                Position pos =  entry.getValue();
                String chromosome = pos.getChromosome().getName();
                int posInt = pos.getPosition();
                RefTagData rtd = new RefTagData(tag,chromosome,posInt,refGenome);
                if(reftagReftagIDMap.containsKey(rtd)) continue;  //it is already in the DB skip                    
                refTagInsertPS.setBytes(1, tag.seq2BitAsBytes());
                refTagInsertPS.setShort(2, tag.seqLength());
                refTagInsertPS.setString(3, pos.getChromosome().getName());
                refTagInsertPS.setInt(4, pos.getPosition());
                refTagInsertPS.setInt(5, refGenomeID);
                
                refTagInsertPS.addBatch();
                batchCount++;
                totalCount++;
                if(batchCount>100000) {
                   // System.out.println("tagInsertPS.executeBatch() "+batchCount);
                    refTagInsertPS.executeBatch();
                    //connection.commit();
                    batchCount=0;
                }
            }
            refTagInsertPS.executeBatch();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
            return false;
        }

        if(totalCount>0) {
            System.out.println("RepGenSQLite:putAllRefTag, totalCount=" + totalCount + ",loadingHash\n");
            loadRefTagHash();
        } 
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

    //This method is called from RepGenAlignerPlugin to add  tags
    // to the  tagAlignments table.
    @Override
    public void putTagTagAlignments(Multimap<Tag,AlignmentInfo> tagAlignInfoMap) {
        int batchCount=0;
        try {
            connection.setAutoCommit(false);
            
            for (Map.Entry<Tag, AlignmentInfo> entry : tagAlignInfoMap.entries()) {
                // Put tag alignments into the tagAlignments table
                AlignmentInfo ai=entry.getValue();
                int ind=1;
 
                tagtagAlignmentInsertPS.setInt(ind++, tagTagIDMap.get(entry.getKey()));
                tagtagAlignmentInsertPS.setInt(ind++, tagTagIDMap.get(ai.tag2()));
                tagtagAlignmentInsertPS.setInt(ind++, ai.score());  // alignment score

                tagtagAlignmentInsertPS.addBatch();
                batchCount++;
                if(batchCount>100000) {
                   // System.out.println("putTagAlignments next"+batchCount);
                    tagtagAlignmentInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            tagtagAlignmentInsertPS.executeBatch();
            connection.setAutoCommit(true);
            // print some metrics for debugging
            ResultSet rs = connection.createStatement().executeQuery("select count (*) as numAlignments from tag_tag_Alignments");
            if (rs.next()) {
                System.out.println("Total alignments in tag_tag_Alignments table: " + rs.getInt("numAlignments"));
            }
 
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    @Override
    public void putTagRefTagAlignments(Multimap<Tag,AlignmentInfo> tagAlignInfoMap, String refGenome) {
        int batchCount=0;
        try {
            connection.setAutoCommit(false);
            for (Map.Entry<Tag, AlignmentInfo> entry : tagAlignInfoMap.entries()) {
                // Put tag alignments into the tagAlignments table
                AlignmentInfo ai=entry.getValue();
                int ind=1;
 
                tagReftagAlignmentInsertPS.setInt(ind++, tagTagIDMap.get(entry.getKey()));
                RefTagData rtd = new RefTagData(ai.tag2(),ai.tag2chrom(),ai.tag2pos(),refGenome);
                tagReftagAlignmentInsertPS.setInt(ind++, reftagReftagIDMap.get(rtd));
                tagReftagAlignmentInsertPS.setInt(ind++, ai.score());  // alignment score
                
                // WHen retrieving data and the ref_strand is 0, user should
                // reverse-complement the refTag sequence to see the alignment.
                tagReftagAlignmentInsertPS.setInt(ind++,ai.alignmentPos()); 
                tagReftagAlignmentInsertPS.setInt(ind++, ai.ref_strand()); 

                tagReftagAlignmentInsertPS.addBatch();
                batchCount++;
                if(batchCount>100000) {
                   // System.out.println("putTagAlignments next"+batchCount);
                    tagReftagAlignmentInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            tagReftagAlignmentInsertPS.executeBatch();
            connection.setAutoCommit(true);
            // print some metrics for debugging
            ResultSet rs = connection.createStatement().executeQuery("select count (*) as numAlignments from tag_reftag_Alignments");
            if (rs.next()) {
                System.out.println("Total alignments in tag_RefTag_Alignments table: " + rs.getInt("numAlignments"));
            }
 
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }
    //This method is called from RepGenAlignerPlugin to add reftag-reftag alignments
    @Override
    public void putRefRefAlignments(Multimap<RefTagData,AlignmentInfo> tagAlignInfoMap, String refGenome) {
        int batchCount=0;
        try {
            connection.setAutoCommit(false);
            for (Map.Entry<RefTagData, AlignmentInfo> entry : tagAlignInfoMap.entries()) {
                // Put tag alignments into the tagAlignments table
                AlignmentInfo ai=entry.getValue();
                int ind=1;
                refRefAlignmentInsertPS.setInt(ind++, reftagReftagIDMap.get(entry.getKey()));
                RefTagData rtd = new RefTagData(ai.tag2(),ai.tag2chrom(),ai.tag2pos(),refGenome);
                refRefAlignmentInsertPS.setInt(ind++, reftagReftagIDMap.get(rtd));
                refRefAlignmentInsertPS.setInt(ind++, ai.score());  // alignment score

                refRefAlignmentInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                   // System.out.println("putTagAlignments next"+batchCount);
                    refRefAlignmentInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            refRefAlignmentInsertPS.executeBatch();
            connection.setAutoCommit(true);
            // print some metrics for debugging
            ResultSet rs = connection.createStatement().executeQuery("select count (*) as numAlignments from reftag_reftag_Alignments");
            if (rs.next()) {
                System.out.println("Total number of alignments in reftag_reftag_alignments table: " + rs.getInt("numAlignments"));
            }
 
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }
    
    //THis method is called from RepGenAlignerPlugin to add reference kmers
    // to the db.  This method adds data to the tag, physicalMapPosition, and
    // tagMapping tables.
    @Override
    public void putRefTagMapping(Multimap<Tag, Position> refTagPositionMap, String refGenome) {
        int batchCount=0;
        loadReferenceGenomeHash();
        try {
            putAllRefTag(refTagPositionMap,refGenome);
            System.out.println("putREfTagMaping: size of map: " + refTagPositionMap.size()
              + ", keyset size: " + refTagPositionMap.keySet().size() + ", values size: " 
                    + refTagPositionMap.values().size());
            putPhysicalMapPositionsIfAbsent(refTagPositionMap.values(),refGenome);
            connection.setAutoCommit(false);
            for (Map.Entry<Tag, Position> entry : refTagPositionMap.entries()) {
                // add reference tags to tagMapping table - map refTag to physicalMapPosition
                Position entrypos=entry.getValue();
                int ind=1;
                RefTagData rtd = new RefTagData(entry.getKey(),entrypos.getChromosome().getName(), 
                        entrypos.getPosition(),refGenome);
               // Defined above:  posTagInsertPS=connection.prepareStatement(
                //"INSERT OR IGNORE into TagMapping (tagid, position_id, method_id, bp_error, cm_error)" +
                // " values(?,?,?,?,?)");
                posTagMappingInsertPS.setInt(ind++, reftagReftagIDMap.get(rtd)); // refTagID
                posTagMappingInsertPS.setInt(ind++, physicalMapPositionToIDMap.get(entrypos)); // position
                posTagMappingInsertPS.setInt(ind++, getMappingApproachID(entrypos)); // method_id
                posTagMappingInsertPS.setInt(ind++, 0);  //todo this needs to be input data (bp_error)
                posTagMappingInsertPS.setFloat(ind++, 0);  //todo this needs to be input data (cm_error)

                posTagMappingInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) { // LCJ - changed from 100000 - writes are REALLY slow when db is big, 79921 tags
                   // System.out.println("putTagAlignments next"+batchCount);
                    posTagMappingInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            posTagMappingInsertPS.executeBatch();
            connection.setAutoCommit(true);
            // print some metrics for debugging
            ResultSet rs = connection.createStatement().executeQuery("select count (DISTINCT physical_position) as numPhysicalSites from physicalMapPosition");
            if (rs.next()) {
                System.out.println("Total number of distinct physical position sites: " + rs.getInt("numPhysicalSites"));
            }
            rs = connection.createStatement().executeQuery("select count (*) as numPhysicalSites from physicalMapPosition");
            if (rs.next()) {
                System.out.println("Total number of physical position sites: " + rs.getInt("numPhysicalSites"));
            }
            PreparedStatement physMapNumFromTCPMP = connection.prepareStatement(
                    "select count(*) as numSites from (select count(*) as tgcnt,physical_position from physicalMapPosition " +
                    "GROUP BY physical_position) where tgcnt=?");
            physMapNumFromTCPMP.setInt(1, 1);// having 1 tag
            rs = physMapNumFromTCPMP.executeQuery();

            if (rs.next()) {
                System.out.println("Number of physical position sites with 1 tag: " + rs.getInt("numSites"));
            }
            physMapNumFromTCPMP.setInt(1, 2);// having 2 tag
            rs = physMapNumFromTCPMP.executeQuery();
            if (rs.next()) {
                System.out.println("Number of physical position sites with 2 tags: " + rs.getInt("numSites"));
            }
            physMapNumFromTCPMP.setInt(1, 3);// having 3 tags
            rs = physMapNumFromTCPMP.executeQuery();
            if (rs.next()) {
                System.out.println("Number of physical position sites with 3 tags: " + rs.getInt("numSites"));
            }

            PreparedStatement cutSiteGreaterThanPS = connection.prepareStatement(
                    "select count(*) as numSites from (select count(*) as tgcnt,physical_position from physicalMapPosition " +
                    "GROUP BY physical_position) where tgcnt>?");
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

                //                    snpQualityInsertPS=connection.prepareStatement(
                //                            "INSERT OR IGNORE into snpQuality (snpid, taxasubset ,avgDepth, minorDepthProp, minor2DepthProp, gapDepthProp, " +
                //                                    "propCovered, propCovered2, taxaCntWithMinorAlleleGE2,minorAlleleFreq, inbredF_DGE2)" +
                //                                    " values(?,?,?,?,?,?,?,?,?,?,?)");
                //                    System.out.println("vals = " + vals.size());
                //                    System.out.println(vals.get("inbredF_DGE2").toString());
                snpQualityInsertPS.setInt(ind++, snpPosToIDMap.get(p));
                snpQualityInsertPS.setString(ind++, taxaSubset);
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("avgDepth",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("minorDepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("minor2DepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("gapDepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("propCovered",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("propCovered2",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("taxaCntWithMinorAlleleGE2",0.0));
                //                    System.out.println("MAF:"+vals.getOrDefault("minorAlleleFreqGE2",-1.0));
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

                //                    snpQualityInsertPS=connection.prepareStatement(
                //                            "INSERT OR IGNORE into snpQuality (snpid, taxasubset ,avgDepth, minorDepthProp, minor2DepthProp, gapDepthProp, " +
                //                                    "propCovered, propCovered2, taxaCntWithMinorAlleleGE2,minorAlleleFreq, inbredF_DGE2)" +
                //                                    " values(?,?,?,?,?,?,?,?,?,?,?)");
                //                    System.out.println("vals = " + vals.size());
                //                    System.out.println(vals.get("inbredF_DGE2").toString());
                snpQualityInsertPS.setInt(ind++, snpPosToIDMap.get(p));
                snpQualityInsertPS.setString(ind++, taxaSubset);
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("avgDepth",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("minorDepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("minor2DepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("gapDepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("propCovered",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("propCovered2",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("taxaCntWithMinorAlleleGE2",0.0));
                //                    System.out.println("MAF:"+vals.getOrDefault("minorAlleleFreqGE2",-1.0));
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

    private int getReferenceGenomeID(String refGenome) throws SQLException{                
        Integer val=referenceGenomeToIDMap.get(refGenome);
        if(val==null) {
            connection.createStatement().executeUpdate("insert into referenceGenome (refname) " +
                    "values('"+refGenome+"')");
            loadReferenceGenomeHash();
            return referenceGenomeToIDMap.get(refGenome);
        } else return val;
    }
    
    @Override
    public void addReferenceGenome(String name) {
        Integer val=mappingApproachToIDMap.get(name);
        if(val==null) {
            try {
                connection.createStatement().executeUpdate("insert or ignore into referenceGenome (refname) " +
                        "values('"+name+"')");
            } catch (SQLException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            loadReferenceGenomeHash();
        } 
    }
    
    @Override
    public void addMappingApproach(String name) {
        Integer val=mappingApproachToIDMap.get(name);
        if(val==null) {
            try {
                connection.createStatement().executeUpdate("insert into mappingApproach (approach, software, protocols) " +
                        "values('"+name+"','unknown','unknown')");
            } catch (SQLException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            loadMappingApproachHash();
        } 
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
            putAllTag(tagAlleleMap.keySet(),null);
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
        allelePairWithTagid1PS.setInt(1,tagTagIDMap.get(tag));
        ResultSet rs=allelePairWithTagid1PS.executeQuery();
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
    public Set<Tag> getTagsForAllele(Position position, byte allele) {
        return getTagsForAllele(new SimpleAllele(allele,position));
    }

    @Override
    public Set<Tag> getTagsForAllele(Allele allele) {
        return null;
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
    public Set<Tag> getTags() {
        return tagTagIDMap.keySet();
    }
    
    @Override
    public Set<RefTagData> getRefTags() {
        return reftagReftagIDMap.keySet();
    }
    
    @Override
    public PositionList getPhysicalMapPositions() {
        if(physicalMapPositionToIDMap == null) loadPhysicalMapPositionHash();
        PositionListBuilder plb=new PositionListBuilder();
        physicalMapPositionToIDMap.keySet().stream()
                .forEach(p -> plb.add(p));
        plb.sortPositions();
        return plb.build();
    }

    @Override
    public PositionList getPhysicalMapPositions(Chromosome chromosome, int firstPosition, int lastPosition) {
        PositionListBuilder plb=new PositionListBuilder();
        plb.addAll(getPositionSubMap(chromosome,firstPosition,lastPosition).keySet());
        return plb.build();
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

//    @Override
//    public PositionList getTagCutPositions(boolean onlyBest) {
//        if(physicalMapPositionToIDMap == null) loadTagMappingHash();
//        PositionListBuilder plb=new PositionListBuilder();
//        physicalMapPositionToIDMap.keySet().stream()
//        //.filter(p -> p.isAnnotatedWithValue("isbest","true"))
//        .forEach(p -> plb.add(p));  //todo only best not implemented here
//        plb.sortPositions();
//        return plb.build();
//    }

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
        if(physicalMapPositionToIDMap==null) loadPhysicalMapPositionHash();
        Position startPos=new GeneralPosition.Builder(chromosome,firstPosition).build();
        if(lastPosition<0) lastPosition=Integer.MAX_VALUE;
        Position lastPos=new GeneralPosition.Builder(chromosome,lastPosition).build();
        return physicalMapPositionToIDMap.subMap(startPos,lastPos);
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

    //        @Override
    //        public Map<Position, Map<Tag, TaxaDistribution>> getCutPositionTagTaxaMapX(Chromosome chromosome, int firstPosition, int lastPosition) {
    //            //consider doing this all by SQL if performance suffers
    //            PositionList pl=getTagCutPositions(chromosome,firstPosition,lastPosition,true);
    //            ImmutableMap.Builder<Position, Map<Tag, TaxaDistribution>> positionMapBuilder=new ImmutableMap.Builder<>();
    //            pl.stream().forEach(p -> positionMapBuilder.put(p,getTagsTaxaMap(p)));
    //            //this is slow as each position is a separate transaction
    //            return positionMapBuilder.build();
    //        }

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

    private void putPhysicalMapPositionsIfAbsent(Collection<Position> positions, String refGenome) {
        try {
            int batchCount=0;
            int positionExists = 0;
            int newPosition = 0;
            if(physicalMapPositionToIDMap==null) loadPhysicalMapPositionHash();
            connection.setAutoCommit(false);
            System.out.println("putPhysicalMapPositionsIfAbsent: size of positions: " + positions.size());
            PreparedStatement posInsertPS=connection.prepareStatement(
                    "INSERT OR IGNORE into physicalMapPosition (reference_genome_id, chromosome, physical_position, strand) values(?,?,?,?)");
            for (Position p : positions) {
                if(physicalMapPositionToIDMap.containsKey(p)) {
                    positionExists++;
                    continue;
                }
                newPosition++;
                posInsertPS.setInt(1, getReferenceGenomeID(refGenome)); 
                posInsertPS.setString(2, p.getChromosome().getName());
                posInsertPS.setInt(3, p.getPosition());
                posInsertPS.setByte(4, p.getStrand());
                posInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                    System.out.println("putPhysicalMapPositionsIfAbsent next"+batchCount);
                    posInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            posInsertPS.executeBatch();
            if(batchCount>0) loadPhysicalMapPositionHash();
            connection.setAutoCommit(true);
            System.out.println("putPhysicalMapPositionsIfAbsent: end, positionExists: " + positionExists
                    + ", newPositions: " + newPosition);
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
        // TODO add refTag clearing when this is revisited fro REpGEn
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
            physicalMapPositionToIDMap = null;
            mappingApproachToIDMap = null;
            loadMappingApproachHash(); // this adds "unknown" to mappingApproachToIDMap
            loadReferenceGenomeHash();
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

    // This method is called  when a user wishes to append data to an existing db from
    // RegGenLoadSeqToDBPlugin;  also called from RepGenLDAnalysis for use in matrix
    // correlation calculations.
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
    public PositionList getTagCutPositions(boolean onlyBest) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public void putAllelePairs(Multimap<Tag, Tuple<Tag, Integer>> tagTagAlignMap) {
        
        int batchCount=0;
        try {
            // Add any tags not currently in the db
            putAllTag(tagTagAlignMap.keySet(),null);
            
            connection.setAutoCommit(false);
            for (Map.Entry<Tag, Tuple<Tag,Integer>> entry : tagTagAlignMap.entries()) {
                Tag tag1 = entry.getKey();
                Tag tag2=entry.getValue().x;
                int score = entry.getValue().y;
                
                int ind=1;
                allelePairInsertPS.setInt(ind++, tagTagIDMap.get(tag1));
                allelePairInsertPS.setInt(ind++, tagTagIDMap.get(tag2));
                allelePairInsertPS.setInt(ind++, score);
                allelePairInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                   // System.out.println("putAllelePairs next"+batchCount);
                    allelePairInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            allelePairInsertPS.executeBatch();
            connection.setAutoCommit(true);
        } catch (Exception exc) {
            exc.printStackTrace();
        }
    }

    @Override
    public Multimap<Allele, Map<Tag, TaxaDistribution>> getAllelesTagTaxaDistForSNP(Position position) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public Multimap<Tag, AlignmentInfo> getTagAlignmentsForTags(List<Tag> tags, int minscore) {
        // This method gets all tag-tag alignments for a list of tags (all data is for non-reference tags)
        ImmutableMultimap.Builder<Tag,AlignmentInfo> tagAIBuilder = ImmutableMultimap.builder();
        loadTagHash(); // get updated tag map
        loadRefTagHash();
 
        try {
            for (Tag tag: tags){
                // For each tag on the list, get all its alignments
                Integer tagID = tagTagIDMap.get(tag);
                if (tagID == null) {
                    // what is best to print out here?
                    // should this return a NULL to indicate an error with the tag?
                    System.out.println("getAlignmentsForTag: no tagID in alignments table for tag: " + tag.sequence());
                    continue;
                }
                // The query specifies that tag1 is not a reference tag
                tagAlignForNonRefTagPS.setInt(1,tagID);
                tagAlignForNonRefTagPS.setInt(2, minscore);
                ResultSet rs = tagAlignForNonRefTagPS.executeQuery();
                while (rs.next()) {
                    int score = rs.getInt("score");
                    Tag tag2 = tagTagIDMap.inverse().get(rs.getInt("tag2id"));
                    AlignmentInfo ai = new AlignmentInfo(tag2,null,-1,-1, -1, null,score);
                    tagAIBuilder.put(tag,ai);
                }              
            }
        } catch (SQLException exc) {
            System.out.println("getTagAlignmentsForTag: caught SQLException attempting to grab alignment data ");
            exc.printStackTrace();
        }
        return tagAIBuilder.build();       
    }

    @Override
    public Multimap<Tag, AlignmentInfo> getAllNonRefTagAlignments(int minscore) {
        // This method returns only the tag-tag alignments for all tags
        ImmutableMultimap.Builder<Tag,AlignmentInfo> tagAIBuilder = ImmutableMultimap.builder();
        loadTagHash(); // get updated tag map
        Set<Tag> tagsToAlign = getTags();
        List<Tag> tagList = new ArrayList<Tag>(tagsToAlign);
        return getTagAlignmentsForTags(tagList,minscore);   // return tag-tag alignments    
    }
    
    @Override
    // This returns just reftag-reftag alignments
    public Multimap<RefTagData, AlignmentInfo> getRefAlignmentsForRefTags(List<RefTagData> refTags, int minscore) {
        ImmutableMultimap.Builder<RefTagData,AlignmentInfo> tagAIBuilder = ImmutableMultimap.builder();
        loadRefTagHash();
        try {
            for (RefTagData tag: refTags){
                // For each tag on the list, get all its alignments
                Integer tagID = reftagReftagIDMap.get(tag);
                if (tagID == null) {
                    // what is best to print out here?
                    System.out.println("getAlignmentsForTag: no tagID in alignments table for reftag: " + tag.tag().sequence());
                    continue;
                }
                // The query specifies that tag1ID IS a ref tag
                refAlignForRefTagPS.setInt(1,tagID);
                refAlignForRefTagPS.setInt(2, minscore);
                ResultSet rs = refAlignForRefTagPS.executeQuery();
                while (rs.next()) {
                    int score = rs.getInt("score");
                    
                    RefTagData tag2data = reftagReftagIDMap.inverse().get(rs.getInt("tag2id"));
                    AlignmentInfo ai = new AlignmentInfo(tag2data.tag(),tag2data.chromosome(),tag2data.position(),-1,1,tag2data.refGenome(),score);
                    tagAIBuilder.put(tag,ai);
                }              
            }
        } catch (SQLException exc) {
            System.out.println("getAlignmentsForRefTags: caught SQLException attempting to grab taxa Distribution ");
            exc.printStackTrace();
        }
        return tagAIBuilder.build();             
    }

    @Override
    public Multimap<RefTagData, AlignmentInfo> getAllRefTagAlignments(int minscore) {
        // for each tag on the refTag list, grab all it's alignments.
        // This is just getting reftag/reftag alignments
        loadRefTagHash(); // get updated list of refTags
 
        Set<RefTagData> refTags = getRefTags();
        List<RefTagData> refTagList = new ArrayList<RefTagData>(refTags);
        return getRefAlignmentsForRefTags(refTagList,minscore);           
    }
    
    @Override
    public Multimap<Tag, AlignmentInfo> getRefAlignmentsForTags(List<Tag> tags, int minscore) {
        // This method gets all ref-tag alignments for each non-ref tag on the list
        ImmutableMultimap.Builder<Tag,AlignmentInfo> tagAIBuilder = ImmutableMultimap.builder();
        // DO we need to reload the hashes?  WIll new tags be added after we've loaded the db?
//        loadTagHash(); // get updated tag map
//        loadRefTagHash();
 
        try {
            for (Tag tag: tags){
                // For each tag on the list, get all its alignments
                Integer tagID = tagTagIDMap.get(tag);
                if (tagID == null) {
                    // what is best to print out here?
                    // should this return a NULL to indicate an error with the tag?
                    System.out.println("getAlignmentsForTag: no tagID in alignments table for tag: " + tag.sequence());
                    continue;
                }
 
                refTagAlignsForNonRefTagPS.setInt(1,tagID);
                refTagAlignsForNonRefTagPS.setInt(2, minscore);
                ResultSet rs = refTagAlignsForNonRefTagPS.executeQuery();
                while (rs.next()) {
                    int alignPos = rs.getInt("ref_align_start_pos");
                    int ref_strand = rs.getInt("ref_align_strand");
                    int score = rs.getInt("score");
                    // refTagID identifies the RefTagData object, that object contains the referenceGenome
                    // name, which is not stored explicitly in the tagAlignments table.
                    RefTagData rtd = reftagReftagIDMap.inverse().get(rs.getInt("refTagID"));
                    AlignmentInfo ai = new AlignmentInfo(rtd.tag(),rtd.chromosome(),
                            rtd.position(),alignPos,ref_strand,rtd.refGenome(),score);                   
                    tagAIBuilder.put(tag,ai);
                }              
            }
        } catch (SQLException exc) {
            System.out.println("getRefAlignmentsForTags: caught SQLException attempting to get reftag alignments ");
            exc.printStackTrace();
        }
        return tagAIBuilder.build();      
    }
    
    @Override
    public Multimap<Tag, AlignmentInfo> getTagAlignmentsForRefTag(RefTagData refTag, int minscore) {
        ImmutableMultimap.Builder<Tag,AlignmentInfo> tagAIBuilder = ImmutableMultimap.builder();
        // This method takes a refTag and finds all the tag (ie non-ref tags) aligned against it
        // that have the specified minimum score
        int refTagID = reftagReftagIDMap.get(refTag);
        try {
            connection.setAutoCommit(false);
            nonReftagAlignmentsForRefTagPS.setInt(1, refTagID );
            nonReftagAlignmentsForRefTagPS.setInt(2, minscore);
            
            ResultSet rs = nonReftagAlignmentsForRefTagPS.executeQuery();
            while (rs.next()) {
                int alignPos = rs.getInt("ref_align_start_pos");
                int ref_strand = rs.getInt("ref_align_strand");
                int score = rs.getInt("score");
                Tag tag1 = tagTagIDMap.inverse().get(rs.getInt("tag1id"));

                AlignmentInfo ai = new AlignmentInfo(refTag.tag(),refTag.chromosome(),
                        refTag.position(),alignPos,ref_strand,refTag.refGenome(),score);                   
                tagAIBuilder.put(tag1,ai);
            }              
        } catch (Exception exc) {
            exc.printStackTrace();
        }
        return tagAIBuilder.build();
    }
    
    @Override
    public void  putTagTagCorrelationMatrix(Multimap<Tag,TagCorrelationInfo> tagCorrelationMap){
        int batchCount=0;
        loadTagHash(); // get updated list of tags
        try {
            connection.setAutoCommit(false);
            for (Map.Entry<Tag, TagCorrelationInfo> entry : tagCorrelationMap.entries()) {
                // Put tag alignments into the tagAlignments table
                TagCorrelationInfo tci=entry.getValue();
                int ind=1;
 
                tagTagCorrelationInsertPS.setInt(ind++, tagTagIDMap.get(entry.getKey()));
                tagTagCorrelationInsertPS.setInt(ind++, tagTagIDMap.get(tci.tag2()));
                tagTagCorrelationInsertPS.setDouble(ind++, tci.t1t2_pearson());
                tagTagCorrelationInsertPS.setDouble(ind++, tci.t1t2_spearman()); 
                tagTagCorrelationInsertPS.setDouble(ind++, tci.pres_abs_pearson());  // presence/absence vector matrix Pearson result
                tagTagCorrelationInsertPS.setDouble(ind++, tci.r2()); // presence/absence vector matrix r-squared results
                
                tagTagCorrelationInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                   // System.out.println("putTagAlignments next"+batchCount);
                    tagTagCorrelationInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            tagTagCorrelationInsertPS.executeBatch();
            connection.setAutoCommit(true);
            // print some metrics for debugging
            ResultSet rs = connection.createStatement().executeQuery("select count (*) as numCorrelations from tagCorrelations");
            if (rs.next()) {
                System.out.println("Total tag-tag correlations in tagCorrelations table: " + rs.getInt("numCorrelations"));
            }
 
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    // A tag may show up as tag1 or tag2 in the tagCorrelations table.
    // for a tagX/tagY correlation, the table will contain either an entry for tag1=x, tag2=Y, or
    // an entry for tag1=y/tag2=x.  The info is the same - only 1 entry is present.
    @Override
    public Multimap<Tag, TagCorrelationInfo> getCorrelationsForTags(List<Tag> tags) {
        ImmutableMultimap.Builder<Tag,TagCorrelationInfo> tagCorBuilder = ImmutableMultimap.builder();
        
        //loadTagHash(); // get updated tag map
        //loadRefTagHash();
 
        try {
//            ResultSet rs = connection.createStatement().executeQuery("select count(*) from tagCorrelations");
//            int size=rs.getInt(1);
//            System.out.println("size of all tagCorrelations in tagCorrelations table=" + size);
            for (Tag tag: tags){
                // For each tag on the list, get all its alignments
                Integer tagID = tagTagIDMap.get(tag);
                if (tagID == null) {
                    // what is best to print out here?
                    // should this return a NULL to indicate an error with the tag?
                    System.out.println("getCorrelationsForTag: no tagID in tagCorrelations table for tag: " + tag.sequence());
                    continue;
                }
                // The query specifies that tag1 is not a reference tag
                tagCorrelationsForTag1PS.setInt(1,tagID);
                ResultSet rs = tagCorrelationsForTag1PS.executeQuery();
                int numCorrelations = 0;
                while (rs.next()) {
                    numCorrelations++;                  
                    int tag2id = rs.getInt("tag2id");
                    double t1t2_p = rs.getFloat("t1t2_pearson");
                    double t1t2_s = rs.getFloat("t1t2_spearman");
                    double pa_pearson = rs.getFloat("pres_abs_pearson");
                    double r2 = rs.getFloat("r2");
                    
                    Tag tag2 = tagTagIDMap.inverse().get(tag2id);
                    TagCorrelationInfo tci = new TagCorrelationInfo(tag2, t1t2_p, t1t2_s, pa_pearson,r2);
 
                    tagCorBuilder.put(tag,tci);
                } 
                //System.out.println("LCJ - RepGenSQLite:getCorrelationsForTag - num correlations found for tag as tag1" + numCorrelations);
                // grab the correlations when tag2 matches the tagID
                tagCorrelationsForTag2PS.setInt(1,tagID);
                rs = tagCorrelationsForTag2PS.executeQuery();
                int numCorrelations2 = 0;
                while (rs.next()) {
                    numCorrelations2++;
                    
                    int tag1id = rs.getInt("tag1id");
                    double t1t2_p = rs.getFloat("t1t2_pearson");
                    double t1t2_s = rs.getFloat("t1t2_spearman");
                    double pa_pearson = rs.getFloat("pres_abs_pearson");
                    double r2 = rs.getFloat("r2");
                    
                    Tag tag2 = tagTagIDMap.inverse().get(tag1id);
                    TagCorrelationInfo tci = new TagCorrelationInfo(tag2, t1t2_p, t1t2_s, pa_pearson,r2);
 
                    tagCorBuilder.put(tag,tci);
                } 
//                System.out.println("LCJ - RepGenSQLite:getCorrelationsForTag - num correlations found for tag as tag1 " + numCorrelations
//                        + ", num correlations found for tag as tag2 " + numCorrelations2);
            }
        } catch (SQLException exc) {
            System.out.println("getAllTaxaMap: caught SQLException attempting to grab taxa Distribution ");
            exc.printStackTrace();
        }
        return tagCorBuilder.build();       
    }
}

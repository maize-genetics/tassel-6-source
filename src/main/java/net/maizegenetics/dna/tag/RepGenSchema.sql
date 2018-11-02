-- Table: tag
CREATE TABLE tag (
    tagid    INTEGER PRIMARY KEY,
    tagName VARCHAR,
    sequence BLOB NOT NULL,
    seqlen INTEGER NOT NULL,
    isReference BOOLEAN,
    qualityScore TEXT,
    numTagInstances INTEGER,
    UNIQUE (sequence, seqlen)
);

-- Table: reftag
CREATE TABLE reftag (
    reftagid    INTEGER PRIMARY KEY,
    tagName VARCHAR,
    sequence BLOB NOT NULL,
    seqlen INTEGER NOT NULL,
    chromosome TEXT,
    position INTEGER,
    refGenomeID INTEGER,
    UNIQUE (sequence, chromosome, seqlen, position, refGenomeID)
);
-- Table: tagMapping
-- Junction (link) table between reftag, and physicalMapPosition
CREATE TABLE tagMapping (
    reftagid       INTEGER NOT NULL,
    position_id  INTEGER NOT NULL,
    method_id    INTEGER NOT NULL,
    bp_error	INTEGER,
    cm_error	FLOAT(2),
    PRIMARY KEY (reftagid, position_id)
);

-- Table: referenceGenome
-- Identifies reference genomes for physicalMapPosition table
CREATE TABLE referenceGenome (
	refid	INTEGER PRIMARY KEY,
	refname	TEXT,
	UNIQUE (refname)
);
 
-- Table: physicalMapPosition
CREATE TABLE physicalMapPosition (
    posid INTEGER   PRIMARY KEY,
    reference_genome_id	INTEGER,
    chromosome TEXT      NOT NULL,
    physical_position   INTEGER   NOT NULL,
    strand     INTEGER(1)  NOT NULL
);
CREATE UNIQUE INDEX physpos_idx ON physicalMapPosition(chromosome,physical_position,strand);
CREATE INDEX phychrpos_idx ON physicalMapPosition(chromosome);

-- Table: tagTagAlignments
CREATE TABLE tag_tag_Alignments (
    tagAlignId INTEGER   PRIMARY KEY,
    tag1id	INTEGER,
    tag2id  INTEGER,
    score INTEGER
);
CREATE INDEX tag_tag_Alignments_idx1 on tag_tag_Alignments(tag1id,tag2id);

-- Table: tagAlignments
-- The ref_align_pos value is where tag1 aligned to tag2 when tag1=non-ref and tag2=ref
-- With non-ref/ref tag comparisons, tag1 is always the non-ref and tag2 is the ref
-- ref_align_position is oftn NOT the start of tag2 as tag2 may have differing
-- length from tag1
CREATE TABLE tag_reftag_Alignments (
    tagAlignId INTEGER   PRIMARY KEY,
    tag1id      INTEGER,
    refTagID  INTEGER,
    score INTEGER,
    ref_align_start_pos  INTEGER,
    ref_align_strand INTEGER
);
CREATE INDEX tag_reftag_Alignments_idx1 on tag_reftag_Alignments(tag1id,score);

-- Table: refTagrefTagAlignments
CREATE TABLE refTag_refTag_Alignments (
    tagAlignId INTEGER   PRIMARY KEY,
    tag1id	INTEGER,
    tag2id  INTEGER,
    score INTEGER
);
CREATE INDEX reftag_reftag_Alignments_idx1 on refTag_refTag_Alignments(tag1id,tag2id);

-- Table: tagCorrelations
-- t1t2_pearson is the pearson correlation of tag1 depths vector by tag2 depths vector
-- t1t2_spearman is the spearman correlation of tag1 depths vector by tag2 depths vector
-- pres_abs_pearson is the pearson correltaion of tag1 prime vec by tag2 prime vec (presence/absence)
-- r2 is the t1' x t2' r2 value from analysis.popgen.LinkageDisequilibrium.calculateRSqr()
CREATE TABLE tagCorrelations (
    tagcorrelationsId INTEGER   PRIMARY KEY,
    tag1id	INTEGER,
    tag2id  INTEGER,
    t1t2_pearson  real,
    t1t2_spearman  real,
    pres_abs_pearson real,
    r2  real
);
CREATE INDEX tagCorrelationsions_idx1 on tagCorrelations(tag1id,tag2id);


-- Table: mappingApproach
CREATE TABLE mappingApproach (
    mapappid INTEGER   PRIMARY KEY,
    approach TEXT NOT NULL UNIQUE,
    software   TEXT NOT NULL,
    protocols   TEXT NOT NULL
);

-- Table: allelepair used for holding Alleles
-- Junction (link) table between tag and alleles
CREATE TABLE allelepair (
    allelepairid        INTEGER PRIMARY KEY,
    tagid_1     INTEGER NOT NULL,
    tagid_2	INTEGER,
    qualityscore INTEGER (1),
    dist_bp	INTEGER,
    dist_cm FLOAT(2)
);
CREATE INDEX newalleleidta_idx on allelepair(tagid_1);

-- Table: tagsnp used for holding Alleles
-- Junction (link) table between tag and snpposition
CREATE TABLE allele (
  alleleid        INTEGER   PRIMARY KEY,
  snpid     INTEGER NOT NULL,
  allelecall         INTEGER(1) NOT NULL,
  qualityscore INTEGER (1)
);
CREATE UNIQUE INDEX snpidallcall_idx on allele (snpid, allelecall);
CREATE INDEX newsnpidallcall_idx on allele (snpid);

-- Table: SNP Position
-- Is the refAllele "NOT NULL" required?
CREATE TABLE snpposition (
    snpid INTEGER   PRIMARY KEY,
    chromosome TEXT      NOT NULL,
    position   INTEGER   NOT NULL,
    strand     INTEGER(1)  NOT NULL,
    qualityscore FLOAT(1),
    refAllele         INTEGER(1) NOT NULL
);
CREATE UNIQUE INDEX snppos_idx ON snpposition(chromosome,position,strand);

-- Table: tagtaxadistribution
CREATE TABLE tagtaxadistribution (
    tagtxdstid  INTEGER   PRIMARY KEY,
    tagid      INTEGER NOT NULL,
    depthsRLE  BLOB,
    totalDepth INTEGER
);
CREATE INDEX tagid_idx ON tagtaxadistribution(tagid);


-- Table: taxa
CREATE TABLE taxa (
    taxonid INTEGER PRIMARY KEY,
    name    TEXT NOT NULL
);

-- Table: SNP Quality
--
CREATE TABLE snpQuality (
  snpid INTEGER   NOT NULL,
  taxasubset TEXT      NOT NULL,
  avgDepth REAL NOT NULL,
  minorDepthProp REAL NOT NULL,
  minor2DepthProp REAL NOT NULL,
  gapDepthProp REAL NOT NULL,
  propCovered REAL NOT NULL,
  propCovered2 REAL NOT NULL,
  taxaCntWithMinorAlleleGE2 REAL NOT NULL,
  minorAlleleFreqGE2      REAL NOT NULL,
  inbredF_DGE2    REAL,
  externalPositive REAL,
  externalNegative REAL
);
CREATE UNIQUE INDEX snpqual_idx ON snpQuality(snpid, taxasubset);


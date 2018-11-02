-- Table: tag
CREATE TABLE tag (
    tagid    INTEGER PRIMARY KEY,
    tagName VARCHAR,
    sequence BLOB NOT NULL,
    seqlen INTEGER NOT NULL,
    UNIQUE (sequence, seqlen)
);

-- Table: tagPosition
-- Junction (link) table between tag, position, and mapping approach
CREATE TABLE tagCutPosition (
    tagid       INTEGER NOT NULL,
    positionid  INTEGER NOT NULL,
    mapappid    INTEGER NOT NULL,
    bestmapping BOOLEAN,
    forward     BOOLEAN,
    cigar       TEXT,
    supportval  INTEGER(2),
    PRIMARY KEY (tagid, positionid, mapappid)
);


-- Table: cutposition
CREATE TABLE cutposition (
    positionid INTEGER   PRIMARY KEY,
    chromosome TEXT      NOT NULL,
    position   INTEGER   NOT NULL,
    strand     INTEGER(1)  NOT NULL
);
CREATE UNIQUE INDEX cutpos_idx ON cutposition(chromosome,position,strand);
CREATE INDEX cutchrpos_idx ON cutposition(chromosome);

-- Table: mappingApproach
CREATE TABLE mappingApproach (
    mapappid INTEGER   PRIMARY KEY,
    approach TEXT NOT NULL UNIQUE,
    software   TEXT NOT NULL,
    protocols   TEXT NOT NULL
);

-- Table: tagsnp used for holding Alleles
-- Junction (link) table between tag and snpposition
CREATE TABLE tagallele (
    tagid        INTEGER NOT NULL,
    alleleid     INTEGER NOT NULL,
    qualityscore INTEGER (1),
    PRIMARY KEY (tagid, alleleid)
);
CREATE INDEX newalleleidta_idx on tagallele(alleleid);

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

-- Table: tissue
CREATE TABLE tissue (
  tissueid INTEGER PRIMARY KEY,
  tissue    TEXT NOT NULL
);

-- Table: tagTaxaTissueDist
CREATE TABLE tagTaxaTissueDist (
  tagtxdstid  INTEGER   PRIMARY KEY,
  tagid      INTEGER NOT NULL,
  tissueid      INTEGER NOT NULL,
  taxonid      INTEGER NOT NULL,
  readCount  INTEGER
);
CREATE INDEX tagtxdstid_idx ON tagtaxatissuedist(tagtxdstid);
CREATE UNIQUE INDEX tagtistaxon_idx ON tagtaxatissuedist(tagid,tissueid,taxonid);

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


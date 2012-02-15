/*=========================PIXELLUMIDATA=======================*/
CREATE TABLE "PIXELLUMIDATA"
(	"DATA_ID" NUMBER(20,0),
	"ENTRY_ID" NUMBER(20,0),
	"ENTRY_NAME" VARCHAR2(56),
	"SOURCE" VARCHAR2(128),
	"RUNNUM" NUMBER(10,0),
	"NOMINALEGEV" BINARY_FLOAT,
	 CONSTRAINT "PIXELLUMIDATA_PK" PRIMARY KEY ("DATA_ID")
)
PARTITION BY RANGE(RUNNUM)
(
PARTITION YEAR2010 VALUES LESS THAN(153368),
PARTITION YEAR2011 VALUES LESS THAN(183338),
PARTITION THEREST VALUES LESS THAN(MAXVALUE)
);
GRANT SELECT ON "PIXELLUMIDATA" TO PUBLIC;
CREATE INDEX PIXELLUMIDATA_RUNNUM_IDX ON PIXELLUMIDATA(RUNNUM) LOCAL;
/*=========================PIXELLUMIDATA_ID=======================*/
CREATE TABLE "PIXELLUMIDATA_ID"
(	"NEXTID" NUMBER(20,0),
	 CONSTRAINT "PIXELLUMIDATA_ID_PK" PRIMARY KEY ("NEXTID")
);
INSERT INTO "PIXELLUMIDATA_ID"("NEXTID") VALUES(0);
GRANT SELECT ON "PIXELLUMIDATA_ID" TO PUBLIC;
/*=========================PIXELLUMIDATA_ENTRIES=======================*/
CREATE TABLE "PIXELLUMIDATA_ENTRIES"
(	"ENTRY_ID" NUMBER(20,0),
	"REVISION_ID" NUMBER(20,0),
	"NAME" VARCHAR2(56),
	 CONSTRAINT "PIXELLUMIDATA_ENTRIES_PK" PRIMARY KEY ("ENTRY_ID")
);
GRANT SELECT ON "PIXELLUMIDATA_ENTRIES" TO PUBLIC;
/*=========================PIXELLUMIDATA_ENTRIES_ID=======================*/
CREATE TABLE "PIXELLUMIDATA_ENTRIES_ID"
(	"NEXTID" NUMBER(20,0),
	 CONSTRAINT "PIXELLUMIDATA_ENTRIES_ID_PK" PRIMARY KEY ("NEXTID")
);
INSERT INTO "PIXELLUMIDATA_ENTRIES_ID"("NEXTID") VALUES(0);
GRANT SELECT ON "PIXELLUMIDATA_ENTRIES_ID" TO PUBLIC;
/*=========================PIXELLUMIDATA_REV=======================*/
CREATE TABLE "PIXELLUMIDATA_REV"
(	"DATA_ID" NUMBER(20,0),
	"REVISION_ID" NUMBER(20,0),
	 CONSTRAINT "PIXELLUMIDATA_REV_PK" PRIMARY KEY("DATA_ID","REVISION_ID")
);
GRANT SELECT ON "PIXELLUMIDATA_REV" TO PUBLIC;
/*=========================PIXELLUMISUMMARYV2=======================*/
CREATE TABLE "PIXELLUMISUMMARYV2"
(	"DATA_ID" NUMBER(20,0),
	"RUNNUM" NUMBER(10,0),
	"LUMILSNUM" NUMBER(10,0),
	"CMSLSNUM" NUMBER(10,0),
	"INSTLUMI" BINARY_FLOAT,
	"INSTLUMIERROR" BINARY_FLOAT,
	"INSTLUMIQUALITY" NUMBER(5,0),
	"BEAMSTATUS" VARCHAR2(28),
	"BEAMENERGY" BINARY_FLOAT,
	"NUMORBIT" NUMBER(10,0),
	"STARTORBIT" NUMBER(10,0),
	 CONSTRAINT "PIXELLUMISUMMARYV2_PK" PRIMARY KEY ("DATA_ID", "LUMILSNUM"))
PARTITION BY RANGE(DATA_ID)(
PARTITION PIXELLUMISUMMARYV2_P1 VALUES LESS THAN(500),
PARTITION PIXELLUMISUMMARYV2_P2 VALUES LESS THAN(1000),
PARTITION PIXELLUMISUMMARYV2_P3 VALUES LESS THAN(1500),
PARTITION PIXELLUMISUMMARYV2_P4 VALUES LESS THAN(2000),
PARTITION THEREST VALUES LESS THAN(MAXVALUE)
);
GRANT SELECT ON "PIXELLUMISUMMARYV2" TO PUBLIC;

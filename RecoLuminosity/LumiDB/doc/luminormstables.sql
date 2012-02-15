/*=========================LUMINORMS=======================*/
begin execute immediate 'drop table LUMINORMS';exception when others then null;end;
CREATE TABLE "LUMINORMS"
("DATA_ID" NUMBER(20,0),
"ENTRY_ID" NUMBER(20,0),
"ENTRY_NAME" VARCHAR2(56),
"SOURCE" VARCHAR2(128),
"RUNNUM" NUMBER(10,0),
"NOMINALEGEV" BINARY_FLOAT,
CONSTRAINT "LUMINORMS_PK" PRIMARY KEY ("DATA_ID")
);
GRANT SELECT ON "LUMINORMS" TO PUBLIC;
/*=========================LUMINORMS_ID=======================*/
begin execute immediate 'drop table LUMINORMS_ID';exception when others then null;end;
CREATE TABLE "LUMINORMS_ID"
(	"NEXTID" NUMBER(20,0),
	 CONSTRAINT "LUMINORMS_ID_PK" PRIMARY KEY ("NEXTID")
);
INSERT INTO "LUMINORMS_ID"("NEXTID") VALUES(0);
GRANT SELECT ON "LUMINORMS_ID" TO PUBLIC;
/*=========================LUMINORMS_ENTRIES=======================*/
begin execute immediate 'drop table LUMINORMS_ENTRIES';exception when others then null;end;
CREATE TABLE "LUMINORMS_ENTRIES"
(	"ENTRY_ID" NUMBER(20,0),
	"REVISION_ID" NUMBER(20,0),
	"NAME" VARCHAR2(56),
	 CONSTRAINT "LUMINORMS_ENTRIES_PK" PRIMARY KEY ("ENTRY_ID")
);
GRANT SELECT ON "LUMINORMS_ENTRIES" TO PUBLIC;
/*=========================LUMINORMS_ENTRIES_ID=======================*/
begin execute immediate 'drop table LUMIDATA_ENTRIES_ID';exception when others then null;end;
CREATE TABLE "LUMINORMS_ENTRIES_ID"
(	"NEXTID" NUMBER(20,0),
	 CONSTRAINT "LUMINORMS_ENTRIES_ID_PK" PRIMARY KEY ("NEXTID")
);
INSERT INTO "LUMINORMS_ENTRIES_ID"("NEXTID") VALUES(0);
GRANT SELECT ON "LUMINORMS_ENTRIES_ID" TO PUBLIC;
/*=========================LUMINORMS_REV=======================*/
begin execute immediate 'drop table LUMINORMS_REV';exception when others then null;end;
CREATE TABLE "LUMINORMS_REV"
(	"DATA_ID" NUMBER(20,0),
	"REVISION_ID" NUMBER(20,0),
	 CONSTRAINT "LUMINORMS_REV_PK" PRIMARY KEY("DATA_ID","REVISION_ID")
);
GRANT SELECT ON "LUMINORMS_REV" TO PUBLIC;

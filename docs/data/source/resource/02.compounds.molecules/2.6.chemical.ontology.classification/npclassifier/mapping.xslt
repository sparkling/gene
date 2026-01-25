<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  NPClassifier -> Chemical Ontology Unified Schema Mapping
  ============================================================
  Source: ./schema.md (NPClassifier Natural Product Classifier)
  Target: ../schema.json (2.6 Chemical Ontology and Classification Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ smiles                      │ canonical_smiles        │ Direct       │
  │ pathway                     │ npc_pathway             │ Direct       │
  │ pathway_score               │ npc_pathway_score       │ Direct       │
  │ superclass                  │ npc_superclass          │ Direct       │
  │ superclass_score            │ npc_superclass_score    │ Direct       │
  │ class                       │ npc_class               │ Direct       │
  │ class_score                 │ npc_class_score         │ Direct       │
  │ isglycoside                 │ is_glycoside            │ Boolean      │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null
  - Scores: undefined → null

  NOTES:
  - NPClassifier predicts biosynthetic pathway (7 types)
  - 35+ superclasses, 200+ classes
  - Confidence scores 0-1 for each level
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    xmlns:local="http://local.functions"
    exclude-result-prefixes="xs fn map array local">

  <xsl:output method="json" indent="yes"/>

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== CHEMICAL STRUCTURE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='smiles'])">
          <fn:null key="canonical_smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="canonical_smiles">
            <xsl:value-of select="fn:string[@key='smiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BIOSYNTHETIC PATHWAY ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='pathway'])">
          <fn:null key="npc_pathway"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="npc_pathway">
            <xsl:value-of select="fn:string[@key='pathway']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='pathway_score']">
          <fn:number key="npc_pathway_score">
            <xsl:value-of select="fn:number[@key='pathway_score']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="npc_pathway_score"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SUPERCLASS ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='superclass'])">
          <fn:null key="npc_superclass"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="npc_superclass">
            <xsl:value-of select="fn:string[@key='superclass']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='superclass_score']">
          <fn:number key="npc_superclass_score">
            <xsl:value-of select="fn:number[@key='superclass_score']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="npc_superclass_score"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CLASS ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='class'])">
          <fn:null key="npc_class"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="npc_class">
            <xsl:value-of select="fn:string[@key='class']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='class_score']">
          <fn:number key="npc_class_score">
            <xsl:value-of select="fn:number[@key='class_score']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="npc_class_score"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GLYCOSIDE FLAG ========== -->

      <xsl:choose>
        <xsl:when test="fn:boolean[@key='isglycoside']">
          <fn:boolean key="is_glycoside">
            <xsl:value-of select="fn:boolean[@key='isglycoside']"/>
          </fn:boolean>
        </xsl:when>
        <xsl:when test="fn:string[@key='isglycoside'] = 'true' or fn:string[@key='isglycoside'] = 'True'">
          <fn:boolean key="is_glycoside">true</fn:boolean>
        </xsl:when>
        <xsl:when test="fn:string[@key='isglycoside'] = 'false' or fn:string[@key='isglycoside'] = 'False'">
          <fn:boolean key="is_glycoside">false</fn:boolean>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="is_glycoside"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">NPClassifier</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='smiles']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>

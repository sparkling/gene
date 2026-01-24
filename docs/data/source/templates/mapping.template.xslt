<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  {{SOURCE_NAME}} -> Unified Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified {{DOMAIN}} Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────┬─────────────────┬──────────────┐
  │ Source Field            │ Target Field    │ Transform    │
  ├─────────────────────────┼─────────────────┼──────────────┤
  │ {{source_id}}           │ id              │ Direct       │
  │ {{source_name}}         │ name            │ Direct       │
  │ {{source_field_1}}      │ {{unified_1}}   │ {{type}}     │
  │ {{source_field_2}}      │ {{unified_2}}   │ {{type}}     │
  │ (computed)              │ {{computed}}    │ Concatenate  │
  └─────────────────────────┴─────────────────┴──────────────┘

  NULL HANDLING:
  - {{field_1}}: "" or "-" → null
  - {{field_2}}: -1 → null

  NOTES:
  - {{Any special transformation notes}}
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

  <!-- ============================================================
       Entry Point: Parse JSON input
       ============================================================ -->
  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <!-- ============================================================
       Main Record Transformation
       ============================================================ -->
  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== DIRECT FIELD MAPPINGS ========== -->

      <!-- Primary identifier -->
      <fn:string key="id">
        <xsl:value-of select="fn:string[@key='{{source_id}}']"/>
      </fn:string>

      <!-- Display name -->
      <fn:string key="name">
        <xsl:value-of select="fn:string[@key='{{source_name}}']"/>
      </fn:string>

      <!-- ========== TRANSFORMED FIELDS ========== -->

      <!-- Example: Remove prefix transformation -->
      <!--
      <fn:string key="chromosome">
        <xsl:value-of select="replace(fn:string[@key='Chromosome'], '^chr', '')"/>
      </fn:string>
      -->

      <!-- Example: Enum mapping -->
      <!--
      <fn:string key="classification">
        <xsl:variable name="source_val" select="fn:string[@key='SourceClassification']"/>
        <xsl:choose>
          <xsl:when test="$source_val = 'ValueA'">unified_a</xsl:when>
          <xsl:when test="$source_val = 'ValueB'">unified_b</xsl:when>
          <xsl:when test="$source_val = 'ValueC'">unified_c</xsl:when>
          <xsl:otherwise>unknown</xsl:otherwise>
        </xsl:choose>
      </fn:string>
      -->

      <!-- ========== NULL HANDLING ========== -->

      <!-- Example: Convert sentinel values to null -->
      <!--
      <xsl:choose>
        <xsl:when test="fn:string[@key='optional_field'] = ('-', '', 'NA')">
          <fn:null key="optional_field"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="optional_field">
            <xsl:value-of select="fn:string[@key='optional_field']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>
      -->

      <!-- Example: Numeric sentinel to null -->
      <!--
      <xsl:choose>
        <xsl:when test="fn:number[@key='numeric_id'] = -1">
          <fn:null key="numeric_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:number key="numeric_id">
            <xsl:value-of select="fn:number[@key='numeric_id']"/>
          </fn:number>
        </xsl:otherwise>
      </xsl:choose>
      -->

      <!-- ========== ARRAY TRANSFORMATIONS ========== -->

      <!-- Example: Transform array of objects -->
      <!--
      <fn:array key="items">
        <xsl:for-each select="fn:array[@key='source_items']/fn:map">
          <fn:map>
            <fn:string key="name">
              <xsl:value-of select="fn:string[@key='item_name']"/>
            </fn:string>
            <fn:string key="id">
              <xsl:value-of select="fn:string[@key='item_id']"/>
            </fn:string>
          </fn:map>
        </xsl:for-each>
      </fn:array>
      -->

      <!-- Example: Split delimited string to array -->
      <!--
      <fn:array key="tags">
        <xsl:for-each select="tokenize(fn:string[@key='tags_csv'], ',')">
          <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
        </xsl:for-each>
      </fn:array>
      -->

      <!-- ========== COMPUTED FIELDS ========== -->

      <!-- Example: Concatenate fields -->
      <!--
      <fn:string key="full_notation">
        <xsl:value-of select="concat(
          fn:string[@key='prefix'], ':',
          fn:string[@key='main'], '-',
          fn:string[@key='suffix']
        )"/>
      </fn:string>
      -->

      <!-- ========== NESTED OBJECT TRANSFORMATION ========== -->

      <!-- Example: Transform nested object -->
      <!--
      <fn:map key="coordinates">
        <fn:string key="chromosome">
          <xsl:value-of select="local:normalize-chromosome(fn:map[@key='location']/fn:string[@key='chr'])"/>
        </fn:string>
        <fn:number key="position">
          <xsl:value-of select="fn:map[@key='location']/fn:number[@key='pos']"/>
        </fn:number>
      </fn:map>
      -->

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">{{SOURCE_NAME}}</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='{{source_id}}']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Reusable Functions
       ============================================================ -->

  <!-- Normalize chromosome: remove "chr" prefix, uppercase -->
  <xsl:function name="local:normalize-chromosome" as="xs:string">
    <xsl:param name="chr" as="xs:string"/>
    <xsl:value-of select="replace(upper-case($chr), '^CHR', '')"/>
  </xsl:function>

  <!-- Convert empty/placeholder to empty string (for null handling) -->
  <xsl:function name="local:nullify" as="xs:string">
    <xsl:param name="value" as="xs:string"/>
    <xsl:choose>
      <xsl:when test="$value = ('-', '.', '', 'NA', 'N/A', 'null')"></xsl:when>
      <xsl:otherwise><xsl:value-of select="$value"/></xsl:otherwise>
    </xsl:choose>
  </xsl:function>

</xsl:stylesheet>

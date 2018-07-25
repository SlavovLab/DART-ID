<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
		xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:loc="http://regis-web.systemsbiology.net/pepXML">
  <xsl:output method="text" />

  <xsl:template match="/*/">
    <xsl:text>Found anything: </xsl:text>
    <xsl:value-of select ="local-name()"/><xsl:text>&#10;</xsl:text>
    <xsl:apply-templates select="loc:msms_pipeline_analysis/loc:analysis_summary" />
  </xsl:template>

  <xsl:template match="loc:msms_pipeline_analysis/loc:analysis_summary">
    <xsl:text>Found analysis_summary</xsl:text>
    <xsl:apply-templates select="loc:peptideprophet_summary" />
  </xsl:template>

  <xsl:template match="loc:peptideprophet_summary">
    <xsl:text>Found peptideprophet_summary</xsl:text>
    <xsl:apply-templates select="loc:mixture_model" />
  </xsl:template>

  <xsl:template match="loc:mixture_model">
    <xsl:text>d </xsl:text>
    <xsl:value-of select="@precursor_ion_charge" /><xsl:text> </xsl:text>
    <xsl:value-of select="@prior_probability" /><xsl:text>&#10;</xsl:text>
  </xsl:template>

</xsl:stylesheet>

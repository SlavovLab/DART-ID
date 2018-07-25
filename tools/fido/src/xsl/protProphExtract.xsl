<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
		xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:loc="http://regis-web.systemsbiology.net/protXML">
  <xsl:output method="text" />

  <xsl:template match="*">
    <xsl:apply-templates select="loc:protein_group" />
  </xsl:template>

  <xsl:template match="loc:protein_group">
<!--    <xsl:text>hello</xsl:text>
-->
    <xsl:value-of select="@probability" />
    <xsl:text> { </xsl:text>
    <xsl:value-of select="loc:protein/@protein_name" />
    <xsl:for-each select="loc:protein/loc:indistinguishable_protein">
      <xsl:text> , </xsl:text><xsl:value-of select="@protein_name" />
    </xsl:for-each>
    <xsl:text> }&#10;</xsl:text>
  </xsl:template>

</xsl:stylesheet>
